function projLAT_patient = projLATs_to_Surface(patient2xCm, patient2yCm, patient2zCm, patient2xSurfaceCm, patient2ySurfaceCm, patient2zSurfaceCm, actTimes2)

    ptCloudSur = pointCloud([patient2xCm, patient2yCm, patient2zCm]);
    
    for i = 1:length(patient2xSurfaceCm)
    
        [indiT,distsT] = findNearestNeighbors(ptCloudSur,[ patient2xSurfaceCm(i), patient2ySurfaceCm(i), patient2zSurfaceCm(i) ], 1);
        projLAT_patient(i) = actTimes2(indiT);
    
    end

end