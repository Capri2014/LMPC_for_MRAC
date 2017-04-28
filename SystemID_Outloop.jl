function SystemID_Outloop(time::Array{Int64,1}, it::Int64, SS::Array{Float64,3}, OldU::Array{Float64,3})

    SSdim_ID = sum(time)
    vector_A1   = zeros(SSdim_ID-it, 3)
    vector_A2   = zeros(SSdim_ID-it, 3)

    vector_b1   = zeros(SSdim_ID-it, 1)
    vector_b2   = zeros(SSdim_ID-it, 1)

    Counter_ID = 1
    for ii = 1:it
        for jj = 1:(time[ii]-1)
            vector_A1[Counter_ID,:] = [ SS[1,jj,ii], SS[2,jj,ii], OldU[1,jj,ii] ]
            vector_b1[Counter_ID,:] = [ SS[1,jj+1,ii] ]
            
            vector_A2[Counter_ID,:] = [ SS[1,jj,ii], SS[2,jj,ii], OldU[1,jj,ii]]
            vector_b2[Counter_ID,:] = [ SS[2,jj+1,ii] ]

            Counter_ID = Counter_ID + 1
        end
    end

    Matrix1 = vector_A1'*vector_A1
    Matrix2 = vector_A2'*vector_A2


    DimentionEye = size(Matrix1,1)
    Row1    = (Matrix1 + 0.00000001*eye(DimentionEye, DimentionEye) ) \ (vector_A1' * vector_b1)
    Row2    = (Matrix2 + 0.00000001*eye(DimentionEye, DimentionEye)) \ (vector_A2' * vector_b2)

    
    MeanEstimate = zeros(2,3)
    MeanEstimate[1,:] = Row1
    MeanEstimate[2,:] = Row2

    MSE1 = 0;
    MSE2 = 0;
    for i = 1:(SSdim_ID-it)
        MSE1 = MSE1 + ( vector_b1[i,:] - *(vector_A1[i,:], Row1) )^2
        MSE2 = MSE2 + ( vector_b2[i,:] - vector_A2[i,1]* Row2[1] )^2
    end
    MSE1 = MSE1^(0.5)/(SSdim_ID-it-1)
    MSE2 = MSE2^(0.5)/(SSdim_ID-it-1)

    MSE = zeros(2)
    MSE[1] = MSE1[1]
    MSE[2] = MSE2[1]

    return MeanEstimate, MSE
    
end