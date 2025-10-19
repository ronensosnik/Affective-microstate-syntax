function edges = emptyEdges()
edges = table([],[],[],[],[],[], 'VariableNames',{'src','tgt','IRR','IRR_CIlo','IRR_CIhi','qFDR'});
edges.srcLabel = strings(0,1); edges.tgtLabel = strings(0,1); edges.Direction = strings(0,1);
end