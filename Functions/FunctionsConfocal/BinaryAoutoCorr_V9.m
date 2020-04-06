function [RCh1 RCh2 RChSum RChCross k]=BinaryAoutoCorr_V9(DataSelectionRules,TTTRCh1Time,TTTRCh2Time,B,ncas,TTTRRes)


%%Autocorrelation 
%[RCh1 Lag]=PointerBasedCorrx_V1(TTTRCh1Time'/TTTRRes,B,ncas,'UnBiased');
[RCh1 Lag]=PointerBasedCorrx_V3(TTTRCh1Time'/TTTRRes,TTTRCh1Time'/TTTRRes,B,ncas,'UnBiased');

%[RCh2 Lag]=PointerBasedCorrx_V1(TTTRCh2Time'/TTTRRes,B,ncas,'UnBiased');
[RCh2 Lag]=PointerBasedCorrx_V3(TTTRCh2Time'/TTTRRes,TTTRCh2Time'/TTTRRes,B,ncas,'UnBiased');

%%Autocorrelation on the Sum
TTTRSumTime=sort([TTTRCh1Time' TTTRCh2Time']);
%[RChSum Lag]=PointerBasedCorrx_V1(TTTRSumTime/TTTRRes,B,ncas,'UnBiased');
[RChSum Lag]=PointerBasedCorrx_V3(TTTRSumTime/TTTRRes,TTTRSumTime/TTTRRes,B,ncas,'UnBiased');

%%CrossCorrelation between chanels
[RChCross Lag]=PointerBasedCorrx_V3(TTTRCh1Time'/TTTRRes,TTTRCh2Time'/TTTRRes,B,ncas,'UnBiased');

k=Lag*TTTRRes;%[s]
%%%%%%%%%%%%%%%%%%%%%

