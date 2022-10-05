function var = RotModule
    TP = types;

    EA  = EulerRot;
    PRV = PRVRot;
    EP  = EPRot;
    CRP = CRPRot;
    MRP = MRPRot;
    DCM = DCMRot;
    
    var.euler2dcm = @EA.euler2dcm;
    var.dcm2euler = @EA.dcm2euler;
    var.euler2om  = @EA.euler2om;
    var.om2euler  = @EA.om2euler;

    var.PRV2dcm   = @PRV.PRV2dcm;
    var.dcm2PRV   = @PRV.dcm2PRV;
    var.PRV2om    = @PRV.PRV2om;
    var.om2PRV    = @PRV.om2PRV;
  
    var.EP2dcm    = @EP.EP2dcm;
    var.dcm2EP    = @EP.dcm2EP;
    var.EP2PRV    = @EP.EP2PRV;
    var.PRV2EP    = @EP.PRV2EP;
    var.SumEP     = @EP.SumEP;
    var.EP2om     = @EP.EP2om;
    var.om2EP     = @EP.om2EP;
    var.EPdiff    = @EP.EPdiff;
  
    var.CRP2dcm   = @CRP.CRP2dcm;
    var.dcm2CRP   = @CRP.dcm2CRP;
    var.CRP2PRV   = @CRP.CRP2PRV;
    var.PRV2CRP   = @CRP.PRV2CRP;
    var.SumCRP    = @CRP.SumCRP;
    var.CRP2om    = @CRP.CRP2om;
    var.om2CRP    = @CRP.om2CRP;
  
    var.MRP2dcm   = @MRP.MRP2dcm;
    var.dcm2MRP   = @MRP.dcm2MRP;
    var.MRP2PRV   = @MRP.MRP2PRV;
    var.PRV2MRP   = @MRP.PRV2MRP;
    var.MRP2CRP   = @MRP.MRP2CRP;
    var.CRP2MRP   = @MRP.CRP2MRP;
    var.SumMRP    = @MRP.SumMRP;
    var.MRP2om    = @MRP.MRP2om;
    var.om2MRP    = @MRP.om2MRP;
    var.MRPdiff   = @MRP.MRPdiff;

    var.DCMmul    = @DCM.DCMmul;
    var.DCMtrs    = @DCM.DCMtrs;
end