
function [SGN, EXP, MNT] = Get_IEEE754(floatIn)
% Extract the IEEE-754 structure from a floating point number
% Inputs: floatIn = complex floating-point input vector
% Output: SGN = complex Sign vector, EXP = complex Exponent vector, 
%         MNT = complex Mantissa vector
% [SGN, EXP, MNT] = Get_IEEE754(floatIn)

    maskExp = 2^8-1;
    shftExp = -23;
    maskMnt = 2^23-1;
    
    X_DIM = size(floatIn,1);
    Y_DIM = size(floatIn,2);

    RL_FLT = real(floatIn);
    SGN    = zeros(X_DIM,Y_DIM,'uint32');
    EXP    = zeros(X_DIM,Y_DIM,'uint32');
    MNT    = zeros(X_DIM,Y_DIM,'uint32');
    for y = 1: Y_DIM
        RL_UINT32 = typecast(RL_FLT(:,y), 'uint32');
        if (~isreal(floatIn))
            IM_FLT = imag(floatIn);
            IM_UINT32 = typecast(IM_FLT(:,y), 'uint32');
            SGN(:,y) = complex(bitget(RL_UINT32, 32),...
                          bitget(IM_UINT32, 32));
            EXP(:,y) = complex(bitand(bitshift(RL_UINT32, shftExp), maskExp),...
                          bitand(bitshift(IM_UINT32, shftExp), maskExp));
            MNT(:,y) = complex(bitand(RL_UINT32, maskMnt),...
                          bitand(IM_UINT32, maskMnt));
        else
            SGN(:,y) = bitget(RL_UINT32, 32);
            EXP(:,y) = bitand(bitshift(RL_UINT32, shftExp), maskExp);
            MNT(:,y) = bitand(RL_UINT32, maskMnt);
        end    
    
    end
    
end
