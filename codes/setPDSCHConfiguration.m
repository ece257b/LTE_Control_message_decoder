function pdsch = setPDSCHConfiguration(enb,dci,RNTI)
    pdsch = [];
    try
        pdsch.RNTI = RNTI;
        pdsch.PRBSet = lteDCIResourceAllocation(enb, dci);
        pdsch.NLayers = enb.CellRefP;   
        pdsch.CSI = 'On';
        pdsch.Modulation = {'QPSK','16QAM','64QAM'};
        pdsch.RV = dci.RV;
        if (enb.CellRefP==1)
            pdsch.TxScheme = 'Port0';
        else
            pdsch.TxScheme = 'TxDiversity';
        end
    catch
        pdsch = [];
        return
    end
end