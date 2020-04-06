function [RCh1 RCh2 RChSum RChCross k]=BinaryAoutoCorr_V11(DataSelectionRules,TTTRCh1Time,TTTRCh2Time,B,ncas,extraTLags,TTTRRes)


%%Autocorrelation 
%[RCh1 Lag]=PointerBasedCorrx_V1(TTTRCh1Time'/TTTRRes,B,ncas,'UnBiased');
[RCh1 Lag]=PointerBasedCorrx_V5(TTTRCh1Time'/TTTRRes,TTTRCh1Time'/TTTRRes,B,ncas,extraTLags,'UnBiased');

%[RCh2 Lag]=PointerBasedCorrx_V1(TTTRCh2Time'/TTTRRes,B,ncas,'UnBiased');
[RCh2 Lag]=PointerBasedCorrx_V5(TTTRCh2Time'/TTTRRes,TTTRCh2Time'/TTTRRes,B,ncas,extraTLags,'UnBiased');

%%Autocorrelation on the Sum
TTTRSumTime=sort([TTTRCh1Time' TTTRCh2Time']);
%[RChSum Lag]=PointerBasedCorrx_V1(TTTRSumTime/TTTRRes,B,ncas,'UnBiased');
[RChSum Lag]=PointerBasedCorrx_V5(TTTRSumTime/TTTRRes,TTTRSumTime/TTTRRes,B,ncas,extraTLags,'UnBiased');

%%CrossCorrelation between chanels
[RChCross Lag]=PointerBasedCorrx_V5(TTTRCh1Time'/TTTRRes,TTTRCh2Time'/TTTRRes,B,ncas,extraTLags,'UnBiased');

k=Lag*TTTRRes;%[s]
%%%%%%%%%%%%%%%%%%%%%

