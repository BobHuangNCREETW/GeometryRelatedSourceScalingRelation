function [pre_S4_Asp,pre_S3_Area]=sub_pre_FinalModel(Mw,dip,Ztor,NZbor,Evtlon,Evtlat)
%%%%%%
%      Source scaling relation submitted to the Seismological Research Letter 2024.
%                         by Bob J.Y. Huang 2024.01.08
%%%%%%
% c0: Constant for M<MBP; c1: slope of M>MBP; 
% MBP: break magnitude for fault rupture reach the bottom of seismogenic zone, which is dip dependent 
c0=0.1139; c1=0.523;
% c2; c3 : intercept and slope for dip dependent MBP  
c2=7.174; c3=-0.011;
% Coefficient for Stage 3, representing depth dependent reletions
c4=0.127;c5=-0.205;
% identify region 
minLON_CA=-125;maxLON_CA=-112;minLAT_CA=31;maxLAT_CA=42;
minLON_TW=120;maxLON_TW=123;minLAT_TW=21.7;maxLAT_TW=25.5;
minLON_JP=128.88;maxLON_JP=143;minLAT_JP=29.79;maxLAT_JP=41.85;
if(Evtlon>=minLON_CA&Evtlon<=maxLON_CA&Evtlat>=minLAT_CA&Evtlat<=maxLAT_CA)
  region='CA';region_flag=1; % 1=CA; 2=TW; 3=JP; 4=other
elseif(Evtlon>=minLON_TW&Evtlon<=maxLON_TW&Evtlat>=minLAT_TW&Evtlat<=maxLAT_TW);
  region='TW';region_flag=2;
elseif(Evtlon>=minLON_JP&Evtlon<=maxLON_JP&Evtlat>=minLAT_JP&Evtlat<=maxLAT_JP);
  region='JP';region_flag=3;
else
  region='other';region_flag=4;
end
% Stage 2 model : A bilinear model with varies break magnitude (dip dependent MBP),
% if M<=MBP, use a constant aspect ratio : 1.3 
% elseif M>MBP, use the second part of the bilinear relation
MBP_ASPS2=c2+(c3)*dip;
if(Mw<=MBP_ASPS2) % Mw smaller than or equal to MBP
  pre_S2_Asp=10^(c0);
elseif(Mw>MBP_ASPS2) % Mw larger than MBP
  pre_S2_Asp=10^(c0+c1.*(Mw-MBP_ASPS2));
end
% predicted aspect ratio in S3 model
pre_S3_Asp=10^(log10(pre_S2_Asp)+c4+(c5*NZbor));
% subregion "others" is for other model using earthquakes out of the three subregions
%RTlarge_CA_S3=0.1746;RTlarge_TW_S3=0.0265;RTlarge_JP_S3=-0.0037;RTlarge_ALL_S3=-0.0363;
RTlarge_CA_S3=0.17;RTlarge_TW_S3=0.03;RTlarge_JP_S3=0;RTlarge_others_S3=-0.04; % use 0 for JP is because the original mean residual is smaller than standard error value
if(region_flag==1) %region_flag: 1=CA; 2=TW; 3=JP; 4=other;
  RT_dip_c6=RTlarge_CA_S3;
elseif(region_flag==2)
  RT_dip_c6=RTlarge_TW_S3;
elseif(region_flag==3)
  RT_dip_c6=RTlarge_JP_S3;
elseif(region_flag==4);
  RT_dip_c6=RTlarge_others_S3;
end
%% applied a smooth function because we did not modeled the RT term of small magnitude earthquakes
smooth_M_range=0.5;
if(Mw<=MBP_ASPS2-smooth_M_range)
  pre_S4_Asp=pre_S3_Asp; %c6=0
elseif(Mw>MBP_ASPS2-smooth_M_range&Mw<MBP_ASPS2+smooth_M_range)
  RT_dip_smooth=(1/(2*smooth_M_range))*RT_dip_c6*(Mw-MBP_ASPS2+smooth_M_range);
  pre_S4_Asp=10.^(log10(pre_S3_Asp)+RT_dip_smooth);
elseif(Mw>=MBP_ASPS2+smooth_M_range)
  pre_S4_Asp=10.^(log10(pre_S3_Asp)+RT_dip_c6);
end
%%%%%
% Area scaling
d0=-0.207;d1=-0.068;d2=1.387;
base_model_logA=Mw-4; % use typical constant 4 as a base model in S1_A(M)
MBP_AreaS2=MBP_ASPS2;
MBPL_AreaS2=MBP_ASPS2+d2;
pre_S2_Area=base_model_logA+max([d0*(max([Mw,MBP_AreaS2])-MBPL_AreaS2)+d1,d1]);
pre_S2_Area=10^pre_S2_Area;
% sur: surface rupture; bur: buried rupture; devided by Ztor = 1 km 
d3_AreaS3_sur=-0.237;d4_AreaS3_sur=0.405;
d3_AreaS3_bur=-0.260;d4_AreaS3_bur=0.295;
if(Ztor<1); % surface rupture fault
  pre_S3_Area=10^(log10(pre_S2_Area)+(d3_AreaS3_sur+d4_AreaS3_sur*NZbor));
elseif(Ztor>=1); % buried fault
  pre_S3_Area=10^(log10(pre_S2_Area)+(d3_AreaS3_bur+d4_AreaS3_bur*NZbor));
end



