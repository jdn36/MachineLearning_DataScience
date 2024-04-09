function [APDduration, indStart, indEnd] = calculateAPD(t,v,threshPerc)
  
   dVdt = diff(v)./diff(t);
   peakV = max(v);
   restV = min(v);
   ampTotal = peakV - restV;
   threshV = peakV - ampTotal*(threshPerc./100);

   indStart = find(dVdt == max(dVdt));
   tStart = t(indStart);

   tOfInterest = t(find(v < threshV ));
   indEnd = (find(tOfInterest > tStart+20));
   tEnd = tOfInterest(indEnd(1));

   indEnd = find(t == tEnd);

   APDduration = tEnd - tStart;

end