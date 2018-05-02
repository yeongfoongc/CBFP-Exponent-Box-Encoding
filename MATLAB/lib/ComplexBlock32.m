classdef ComplexBlock32
    %   ComplexBlock32  Declare a numeric class that performs complex block 
    %   floating-point representation.
    %
    %   Input Parameters: complex floating-point vector in single-precision (32-bit), 
    %   logical flag (false/true) indicates whether common exponent encoding/ 
    %   exponent box encoding is chosen.
    %
    %   obj = ComplexBlock32(A, false) quantizes a complex floating-point 
    %   vector into the new class with common exponent encoding.
    %
    %   obj = ComplexBlock32(A, true) quantizes a complex floating-point 
    %   vector into the new class with exponent box encoding.
    %
    properties (GetAccess = private)
        SharedExponent    % 8-bit exponent applied to the block
        SignVector        % 1-bit complex-valued vectors
        MantissaVector    % 23-bit complex-valued vectors
        ShiftVector       % 1-bit complex-valued vectors (applicable for exponent box encoding)
        EncodeMethod      % 1-bit scalar from input parameter initialization
        
        ExponentBias      % 7-bit value exponent bias
        MantissaBitWidth  % value = 24 (IEEE-754 plus 1 hidden bit)
    end 
    
    properties (Constant, GetAccess = private)
        CommonExponentEncode = false;   % flags for checking input
        ExponentBoxEncode = true;
        SPMantissa = 23;                % single-precision
        SPExponent = 8;
    end
    
    properties
        ref              % reference IEEE-754 complex-valued vector
    end
    
    methods
        
        function obj = ComplexBlock32( val, encodingMethod )
        %     ComplexBlock32 Declare a numeric class that performs complex block
        %     floating-point representation
        %
        %     obj = ComplexBlock32(A, false) quantizes a complex floating-point 
        %     vector into the new class with common exponent encoding.
        %
        %     obj = ComplexBlock32(A, true) quantizes a complex floating-point 
        %     vector into the new class with exponent box encoding.
        %            
            if nargin == 2
                if isa(val,'single') && isscalar(encodingMethod) && islogical(encodingMethod)
                    % convert to row vectors
                    if (size(val,2) < size(val,1))
                        val = val.';
                    end
                    obj.ref = [real(val);imag(val)]';
                    obj.ExponentBias = 127;
                    obj.MantissaBitWidth = 23;
                    [S, E, M] = Get_IEEE754(val);
                    obj.SignVector = logical([real(S);imag(S)]');
                    obj.SharedExponent = int32(max([real(E(:)); imag(E(:))]));
                    % default to common exponent encoding
                    if encodingMethod ~= obj.CommonExponentEncode && encodingMethod ~= obj.ExponentBoxEncode
                        encodingMethod = obj.CommonExponentEncode;
                    end
                    obj.EncodeMethod = encodingMethod;
                    mant = int32([real(M); imag(M)]');
                    E = int32(E);
                    ExponentVector = [real(E(:))'; imag(E(:))']';
                    if (obj.EncodeMethod == obj.ExponentBoxEncode)
                        obj.ShiftVector = ExponentVector + obj.MantissaBitWidth <= obj.SharedExponent;
                        ExponentVector(obj.ShiftVector) = ExponentVector(obj.ShiftVector) + obj.MantissaBitWidth;
                    end
                    expDiff = obj.SharedExponent - ExponentVector;
                    obj.MantissaVector = bitshift(bitset(mant,obj.MantissaBitWidth+1), -expDiff);
                else
                    if (~isa(val,'single'))
                        error('Error: Complex value vector must be single-precision floating-point.');
                    end
                    if (~isscalar(encodingMethod))
                        error('Error: Encoding method must be scalar, Example: Common Exponent Encoding ( 0 ), Exponent Box Encoding ( 1 )');
                    end
                    if (~islogical(encodingMethod))
                        error('Error: Encoding method must be logical, Example: Common Exponent Encoding ( 0 ), Exponent Box Encoding ( 1 )');
                    end
                end
            elseif nargin == 1
                error('Error: Complex input vector not initialized, Example: ComplexBlock32( ComplexVector, EncodingMethod )');
            end
        end
        
        function val = encodedValue(obj)
        %     Display quantized numeric values in the complex block representation
        
            if (obj.EncodeMethod == obj.ExponentBoxEncode)
                obj.MantissaVector(obj.ShiftVector) = bitshift(obj.MantissaVector(obj.ShiftVector),-obj.MantissaBitWidth);
            end
            val = (-1).^(single(obj.SignVector)).*(single(obj.MantissaVector)/2^obj.MantissaBitWidth).*2.^(single(obj.SharedExponent-obj.ExponentBias));
            val = complex(val(:,1), val(:,2)).';
        end
        
        function val = encodedMantissaVector(obj)
        %     Display quantized 24-bit mantissa bits vector in complex block representation
            val = obj.MantissaVector;
        end
        
        function val = encodedShiftVector(obj)
        %     Display 1-bit exponent box bits vector in complex block representation
            val = obj.ShiftVector;
        end
        
        function err = encodingError(obj)
        %     Display quantization error between complex block
        %     representation and reference IEEE-754 vector
            err = obj.ref - obj.encodedValue;
        end
                
        function obj3 = plus(obj1, obj2)
        %     + Plus.
        %     X + Y adds matrices X and Y. X and Y must have compatible sizes. In the
        %     simplest cases, they can be the same size or one can be a scalar. Two
        %     inputs have compatible sizes if, for every dimension, the dimension
        %     sizes of the inputs are either the same or one of them is 1.
        %  
        %     C = plus(A,B) is called for the syntax 'A + B' when A or B is an
        %     object.
            OutputExponent = max(obj1.SharedExponent, obj2.SharedExponent);
            MantVec1 = int32(obj1.MantissaVector);
            MantVec2 = int32(obj2.MantissaVector);
            ExponentDiff1 = int32(OutputExponent - obj1.SharedExponent);
            ExponentDiff2 = int32(OutputExponent - obj2.SharedExponent);
            MantVec1 = bitshift(MantVec1, -ExponentDiff1);
            MantVec2 = bitshift(MantVec2, -ExponentDiff2);
            
            if ( obj1.EncodeMethod == obj1.ExponentBoxEncode )
                MantVec1(obj1.ShiftVector) = bitshift(MantVec1(obj1.ShiftVector), -obj1.MantissaBitWidth);
            end
            if ( obj2.EncodeMethod == obj2.ExponentBoxEncode )
                MantVec2(obj2.ShiftVector) = bitshift(MantVec2(obj2.ShiftVector), -obj2.MantissaBitWidth);
            end
            
            MantVec1(obj1.SignVector) = -MantVec1(obj1.SignVector);
            MantVec2(obj2.SignVector) = -MantVec2(obj2.SignVector);
            OutputMantissa = MantVec1 + MantVec2;
            OutputSign = OutputMantissa < 0;
            OutputMantissa = uint32(abs(OutputMantissa));
            
            WordLengthExpand = int32(floor(log2(single(OutputMantissa))) - obj1.MantissaBitWidth);
            OutputMantissa = bitshift( int32(OutputMantissa), -max(WordLengthExpand(:)) );
            OutputExponent = uint32(int32(OutputExponent) + max(WordLengthExpand(:)));
            
            res = (-1).^single(OutputSign).*(single(OutputMantissa)/2^obj1.MantissaBitWidth).*2^(single(OutputExponent)-obj1.ExponentBias);
            obj3 = ComplexBlock32(complex(res(:,1),res(:,2)),obj1.EncodeMethod);
        end
        
        function obj3 = times(obj1, obj2)
        % .*  Array multiply.
        %     X.*Y denotes element-by-element multiplication. X and Y must have
        %     compatible sizes. In the simplest cases, they can be the same size or
        %     one can be a scalar. Two inputs have compatible sizes if, for every
        %     dimension, the dimension sizes of the inputs are either the same or one
        %     of them is 1.
        %  
        %     C = times(A,B) is called for the syntax 'A .* B' when A or B is an
        %     object.    
            sgn_r = [obj1.SignVector(:,1) ; obj2.SignVector(:,1)];
            sgn_i = [obj1.SignVector(:,2) ; obj2.SignVector(:,2)];
            
            if (obj1.EncodeMethod == obj1.ExponentBoxEncode && obj2.EncodeMethod == obj2.ExponentBoxEncode)
                box_r = [obj1.ShiftVector(:,1) ; obj2.ShiftVector(:,1)];
                box_i = [obj1.ShiftVector(:,2) ; obj2.ShiftVector(:,2)];
                
                shft_man_rA = box_r(1:length(box_r)/2) + box_r(length(box_r)/2+1:end);
                shft_man_rB = box_i(1:length(box_i)/2) + box_i(length(box_i)/2+1:end);
                
                shft_man_r = min( shft_man_rA, shft_man_rB );
                
                shft_man_rA = (shft_man_r - shft_man_rA)*obj1.MantissaBitWidth;
                shft_man_rB = (shft_man_r - shft_man_rB)*obj1.MantissaBitWidth;
                
                shft_man_iA = box_r(1:length(box_r)/2) + box_i(length(box_i)/2+1:end);
                shft_man_iB = box_i(1:length(box_i)/2) + box_r(length(box_r)/2+1:end);
                
                shft_man_i = min( shft_man_iA, shft_man_iB );
                
                shft_man_iA = (shft_man_i - shft_man_iA)*obj1.MantissaBitWidth;
                shft_man_iB = (shft_man_i - shft_man_iB)*obj1.MantissaBitWidth;
            end

            man_r = [obj1.MantissaVector(:,1) ; obj2.MantissaVector(:,1)];
            man_i = [obj1.MantissaVector(:,2) ; obj2.MantissaVector(:,2)];

            e_int = int32( obj1.SharedExponent + obj2.SharedExponent - obj1.ExponentBias );
            
            res_man_rA = uint64( man_r(1:length(man_r)/2) ) .* uint64( man_r(length(man_r)/2+1:end) );
            s_rA = sgn_r(1:length(sgn_r)/2) ~= sgn_r(length(sgn_r)/2+1:end);
            
            res_man_rB = uint64( man_i(1:length(man_i)/2) ) .* uint64( man_i(length(man_i)/2+1:end) );
            s_rB = sgn_i(1:length(sgn_i)/2) ~= sgn_i(length(sgn_i)/2+1:end);
                        
            res_man_iA = uint64( man_r(1:length(man_r)/2) ) .* uint64( man_i(length(man_i)/2+1:end) );
            s_iA = sgn_r(1:length(sgn_r)/2) ~= sgn_i(length(sgn_i)/2+1:end);
            
            res_man_iB = uint64( man_i(1:length(man_i)/2) ) .* uint64( man_r(length(man_r)/2+1:end) );
            s_iB = sgn_i(1:length(sgn_i)/2) ~= sgn_r(length(sgn_r)/2+1:end);

            if (obj1.EncodeMethod == obj1.ExponentBoxEncode && obj2.EncodeMethod == obj2.ExponentBoxEncode)
                e_int_r = e_int - int32(shft_man_r).*obj1.MantissaBitWidth;
                res_man_rA = bitshift( res_man_rA, shft_man_rA );
                res_man_rB = bitshift( res_man_rB, shft_man_rB );
                e_int_i = e_int - int32(shft_man_i).*obj1.MantissaBitWidth;
                res_man_iA = bitshift( res_man_iA, shft_man_iA );
                res_man_iB = bitshift( res_man_iB, shft_man_iB );
            end
            
            res_man_r = (-1).^int64(s_rA).*int64(res_man_rA) - (-1).^int64(s_rB).*int64(res_man_rB);
            res_man_i = (-1).^int64(s_iA).*int64(res_man_iA) + (-1).^int64(s_iB).*int64(res_man_iB);
            
            sgn_res_r = res_man_r < 0;
            res_man_r(sgn_res_r) = -res_man_r(sgn_res_r);
            sgn_res_i = res_man_i < 0;
            res_man_i(sgn_res_i) = -res_man_i(sgn_res_i);
            
            if ( obj1.EncodeMethod == obj1.ExponentBoxEncode && obj2.EncodeMethod == obj2.ExponentBoxEncode )
                e3 = max( max( e_int_r ), max( e_int_i ) );
                res_man_r = bitshift( res_man_r, - e3 + e_int_r );
                res_man_i = bitshift( res_man_i, - e3 + e_int_i );
            end
            max_val = max( max( abs(res_man_r) ), max( abs(res_man_i) ) );
            wordlengthExpand = int32(log2(single(max_val))) - obj1.MantissaBitWidth*2;
            shft_all = int32(log2(single(max_val)) - (obj1.MantissaBitWidth)) - wordlengthExpand;

            res_man_r = bitshift( res_man_r, -int32(shft_all) );
            res_man_i = bitshift( res_man_i, -int32(shft_all) );

            e3 = e_int;

            res = complex( (-1).^sgn_res_r.*(single(res_man_r)/2^(obj1.MantissaBitWidth)).*2.^(single(e3) - obj1.ExponentBias), ...
                           (-1).^sgn_res_i.*(single(res_man_i)/2^(obj1.MantissaBitWidth)).*2.^(single(e3) - obj1.ExponentBias) ...
                          );
              
            obj3 = ComplexBlock32( res, obj1.EncodeMethod );
        end
        
        function obj3 = mtimes(obj1, obj2)
        % *   Matrix multiply.
        %     X*Y is the matrix product of X and Y.  Any scalar (a 1-by-1 matrix)
        %     may multiply anything.  Otherwise, the number of columns of X must
        %     equal the number of rows of Y.
        %  
        %     C = mtimes(A,B) is called for the syntax 'A * B' when A or B is an
        %     object.                
            sgn_r = [obj1.SignVector(:,1) ; obj2.SignVector(:,1)];
            sgn_i = [obj1.SignVector(:,2) ; obj2.SignVector(:,2)];
            
            if (obj1.EncodeMethod == obj1.ExponentBoxEncode && obj2.EncodeMethod == obj2.ExponentBoxEncode)
                box_r = [obj1.ShiftVector(:,1) ; obj2.ShiftVector(:,1)];
                box_i = [obj1.ShiftVector(:,2) ; obj2.ShiftVector(:,2)];
                
                shft_man_rA = box_r(1:length(box_r)/2) + box_r(length(box_r)/2+1:end);
                shft_man_rB = box_i(1:length(box_i)/2) + box_i(length(box_i)/2+1:end);
                
                shft_man_r = min( shft_man_rA, shft_man_rB );
                
                shft_man_rA = (shft_man_r - shft_man_rA)*obj1.MantissaBitWidth;
                shft_man_rB = (shft_man_r - shft_man_rB)*obj1.MantissaBitWidth;
                
                shft_man_iA = box_r(1:length(box_r)/2) + box_i(length(box_i)/2+1:end);
                shft_man_iB = box_i(1:length(box_i)/2) + box_r(length(box_r)/2+1:end);
                
                shft_man_i = min( shft_man_iA, shft_man_iB );
                
                shft_man_iA = (shft_man_i - shft_man_iA)*obj1.MantissaBitWidth;
                shft_man_iB = (shft_man_i - shft_man_iB)*obj1.MantissaBitWidth;
            end

            man_r = [obj1.MantissaVector(:,1) ; obj2.MantissaVector(:,1)];
            man_i = [obj1.MantissaVector(:,2) ; obj2.MantissaVector(:,2)];

            e_int = int32( obj1.SharedExponent + obj2.SharedExponent - obj1.ExponentBias );
            
            res_man_rA = uint64( man_r(1:length(man_r)/2) ) .* uint64( man_r(length(man_r)/2+1:end) );
            s_rA = sgn_r(1:length(sgn_r)/2) ~= sgn_r(length(sgn_r)/2+1:end);
            
            res_man_rB = uint64( man_i(1:length(man_i)/2) ) .* uint64( man_i(length(man_i)/2+1:end) );
            s_rB = sgn_i(1:length(sgn_i)/2) ~= sgn_i(length(sgn_i)/2+1:end);
                        
            res_man_iA = uint64( man_r(1:length(man_r)/2) ) .* uint64( man_i(length(man_i)/2+1:end) );
            s_iA = sgn_r(1:length(sgn_r)/2) ~= sgn_i(length(sgn_i)/2+1:end);
            
            res_man_iB = uint64( man_i(1:length(man_i)/2) ) .* uint64( man_r(length(man_r)/2+1:end) );
            s_iB = sgn_i(1:length(sgn_i)/2) ~= sgn_r(length(sgn_r)/2+1:end);

            if (obj1.EncodeMethod == obj1.ExponentBoxEncode && obj2.EncodeMethod == obj2.ExponentBoxEncode)
                e_int_r = e_int - int32(shft_man_r).*obj1.MantissaBitWidth;
                res_man_rA = bitshift( res_man_rA, shft_man_rA );
                res_man_rB = bitshift( res_man_rB, shft_man_rB );
                e_int_i = e_int - int32(shft_man_i).*obj1.MantissaBitWidth;
                res_man_iA = bitshift( res_man_iA, shft_man_iA );
                res_man_iB = bitshift( res_man_iB, shft_man_iB );
            end
            
            res_man_r = (-1).^int64(s_rA).*int64(res_man_rA) - (-1).^int64(s_rB).*int64(res_man_rB);
            res_man_i = (-1).^int64(s_iA).*int64(res_man_iA) + (-1).^int64(s_iB).*int64(res_man_iB);
            
            if ( obj1.EncodeMethod == obj1.ExponentBoxEncode && obj2.EncodeMethod == obj2.ExponentBoxEncode )
                e3 = max( max( e_int_r ), max( e_int_i ) );
                res_man_r = bitshift( res_man_r, - e3 + ( e_int_r ) );
                res_man_i = bitshift( res_man_i, - e3 + ( e_int_i ) );
            end
            
            res_man_r = sum ( res_man_r );
            res_man_i = sum ( res_man_i );

            sgn_res_r = res_man_r < 0;
            res_man_r(sgn_res_r) = -res_man_r(sgn_res_r);
            sgn_res_i = res_man_i < 0;
            res_man_i(sgn_res_i) = -res_man_i(sgn_res_i);
            
            max_val = max( max( abs(res_man_r) ), max( abs(res_man_i) ) );
            wordlengthExpand = int32(log2(single(max_val))) - obj1.MantissaBitWidth*2;
            shft_all = int32(log2(single(max_val)) - (obj1.MantissaBitWidth)) - wordlengthExpand;

            res_man_r = bitshift( res_man_r, -int32(shft_all) );
            res_man_i = bitshift( res_man_i, -int32(shft_all) );

            e3 = e_int;
            
            res = complex( (-1).^sgn_res_r.*(single(res_man_r)/2^(obj1.MantissaBitWidth)).*2.^(single(e3) - obj1.ExponentBias), ...
                           (-1).^sgn_res_i.*(single(res_man_i)/2^(obj1.MantissaBitWidth)).*2.^(single(e3) - obj1.ExponentBias) ...
                          );
              
            obj3 = ComplexBlock32( res, obj1.EncodeMethod );
        end
        
        function obj2 = sum_fft_terms( obj, N, r, k )
        %  sum_fft_terms. Inner product of fft complex exponentials and complex input vector.
        %     S = sum_fft_terms( obj, N, r, k ) is the inner product of the obj and fft complex exponentials. 
        %     obj is an object of class ComplexBlock32
        %     N is scalar value that specifies the vector length
        %     k is index term of the corresponding fft output coefficient
        
            if ( ~isscalar(k) || ~isscalar(r) || ~(mod(k,1)==0) || ~(mod(r,1)==0) )
                error('Error: r and k must be scalar and integer, r: radix-number, k: integer');
            end

            n = single( 0:N-1 ).';

            exponentials = exp(-1j/N*2*pi*n*k);
            
            Exp_Block = ComplexBlock32( exponentials, obj.EncodeMethod );

            obj2 = Exp_Block * obj;
                        
        end
        
        function obj2 = fft(obj, N)
        % fft Discrete Fourier transform.
        %     fft(X) is the discrete Fourier transform (DFT) of vector X.  For
        %     matrices, the fft operation is applied to each column. For N-D
        %     arrays, the fft operation operates on the first non-singleton
        %     dimension.
        %  
        %     fft(X,N) is the N-point fft, padded with zeros if X has less
        %     than N points and truncated if it has more.
        %  
        %     For length N input vector x, the DFT is a length N vector X,
        %     with elements
        %                      N
        %        X(k) =       sum  x(n)*exp(-j*2*pi*(k-1)*(n-1)/N), 1 <= k <= N.
        %                     n=1
        %     The inverse DFT (computed by IFFT) is given by
        %                      N
        %        x(n) = (1/N) sum  X(k)*exp( j*2*pi*(k-1)*(n-1)/N), 1 <= n <= N.
        %                     k=1
            if (nargin == 2)
                if (~(N > 0) || ~isscalar(N) || ~(mod(log2(N)/log2(4),1)==0) )
                    error('Error: N has to be greater than 0, scalar, and base-4 number.');
                end
            elseif (nargin == 1)
                N = 4^ceil( log2( size(obj.ref,1) )/ log2(4) );
            else
                error('fft supports 1-input: complex vector and 2-input: complex vector, N-point FFT');
            end
            if ( N == 4 )
                
                temp1 = sum_fft_terms( obj, N, N, 0 );
                temp2 = sum_fft_terms( obj, N, N, 1 );
                temp3 = sum_fft_terms( obj, N, N, 2 );
                temp4 = sum_fft_terms( obj, N, N, 3 );
                
                intermed1 = [temp1.encodedValue, ...
                             temp2.encodedValue, ...
                             temp3.encodedValue, ...
                             temp4.encodedValue];
                                
                obj2 = ComplexBlock32( intermed1, obj.EncodeMethod );
                
            elseif ( N == 16 )
                
                tempA = obj.encodedValue;
                
                tempA1 = fft( ComplexBlock32( tempA(1:4:13), obj.EncodeMethod ) );
                tempA2 = fft( ComplexBlock32( tempA(2:4:14), obj.EncodeMethod ) );
                tempA3 = fft( ComplexBlock32( tempA(3:4:15), obj.EncodeMethod ) );
                tempA4 = fft( ComplexBlock32( tempA(4:4:16), obj.EncodeMethod ) );
                
                tempB = zeros( size( tempA ) );
                
                tempB(1:4:13) = tempA1.encodedValue;
                tempB(2:4:14) = tempA2.encodedValue;
                tempB(3:4:15) = tempA3.encodedValue;
                tempB(4:4:16) = tempA4.encodedValue;
                
                k = single( (0:3) );
                tempB1 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*0 ) .* tempB(1:4), obj.EncodeMethod ) );
                tempB2 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*1 ) .* tempB(5:8), obj.EncodeMethod ) );
                tempB3 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*2 ) .* tempB(9:12), obj.EncodeMethod ) );
                tempB4 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*3 ) .* tempB(13:16), obj.EncodeMethod ) );
                
                tempC = zeros( size( tempB ) );
                tempC(1:4:13) = tempB1.encodedValue;
                tempC(2:4:14) = tempB2.encodedValue;
                tempC(3:4:15) = tempB3.encodedValue;
                tempC(4:4:16) = tempB4.encodedValue;
                
                obj2 = ComplexBlock32( single( tempC ), obj.EncodeMethod );                
            elseif ( N == 64 )
                
                tempA = obj.encodedValue;
                
                tempA1 = fft( ComplexBlock32( tempA(1:4:61), obj.EncodeMethod ) );
                tempA2 = fft( ComplexBlock32( tempA(2:4:62), obj.EncodeMethod ) );
                tempA3 = fft( ComplexBlock32( tempA(3:4:63), obj.EncodeMethod ) );
                tempA4 = fft( ComplexBlock32( tempA(4:4:64), obj.EncodeMethod ) );
                
                tempB = zeros( size( tempA ) );
                
                tempB(1:4:61) = tempA1.encodedValue;
                tempB(2:4:62) = tempA2.encodedValue;
                tempB(3:4:63) = tempA3.encodedValue;
                tempB(4:4:64) = tempA4.encodedValue;
                
                k = single( (0:3) );
                tempB1 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*0 ) .* tempB(1:4), obj.EncodeMethod ) );
                tempB2 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*1 ) .* tempB(5:8), obj.EncodeMethod ) );
                tempB3 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*2 ) .* tempB(9:12), obj.EncodeMethod ) );
                tempB4 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*3 ) .* tempB(13:16), obj.EncodeMethod ) );
                tempB5 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*4 ) .* tempB(17:20), obj.EncodeMethod ) );
                tempB6 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*5 ) .* tempB(21:24), obj.EncodeMethod ) );
                tempB7 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*6 ) .* tempB(25:28), obj.EncodeMethod ) );
                tempB8 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*7 ) .* tempB(29:32), obj.EncodeMethod ) );
                tempB9 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*8 ) .* tempB(33:36), obj.EncodeMethod ) );
                tempB10 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*9 ) .* tempB(37:40), obj.EncodeMethod ) );
                tempB11 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*10 ) .* tempB(41:44), obj.EncodeMethod ) );
                tempB12 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*11 ) .* tempB(45:48), obj.EncodeMethod ) );
                tempB13 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*12 ) .* tempB(49:52), obj.EncodeMethod ) );
                tempB14 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*13 ) .* tempB(53:56), obj.EncodeMethod ) );
                tempB15 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*14 ) .* tempB(57:60), obj.EncodeMethod ) );
                tempB16 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*15 ) .* tempB(61:64), obj.EncodeMethod ) );
                
                tempC = zeros( size( tempB ) );
                tempC(1:16:49) = tempB1.encodedValue;
                tempC(2:16:50) = tempB2.encodedValue;
                tempC(3:16:51) = tempB3.encodedValue;
                tempC(4:16:52) = tempB4.encodedValue;
                tempC(5:16:53) = tempB5.encodedValue;
                tempC(6:16:54) = tempB6.encodedValue;
                tempC(7:16:55) = tempB7.encodedValue;
                tempC(8:16:56) = tempB8.encodedValue;
                tempC(9:16:57) = tempB9.encodedValue;
                tempC(10:16:58) = tempB10.encodedValue;
                tempC(11:16:59) = tempB11.encodedValue;
                tempC(12:16:60) = tempB12.encodedValue;
                tempC(13:16:61) = tempB13.encodedValue;
                tempC(14:16:62) = tempB14.encodedValue;
                tempC(15:16:63) = tempB15.encodedValue;
                tempC(16:16:64) = tempB16.encodedValue;
                
                obj2 = ComplexBlock32( single( tempC ), obj.EncodeMethod );                
            elseif ( N == 256 )
                
                tempA = obj.encodedValue;
                
                tempA1 = fft( ComplexBlock32( tempA(1:4:253), obj.EncodeMethod ) );
                tempA2 = fft( ComplexBlock32( tempA(2:4:254), obj.EncodeMethod ) );
                tempA3 = fft( ComplexBlock32( tempA(3:4:255), obj.EncodeMethod ) );
                tempA4 = fft( ComplexBlock32( tempA(4:4:256), obj.EncodeMethod ) );
                
                tempB = zeros( size( tempA ) );
                
                tempB(1:4:253) = tempA1.encodedValue;
                tempB(2:4:254) = tempA2.encodedValue;
                tempB(3:4:255) = tempA3.encodedValue;
                tempB(4:4:256) = tempA4.encodedValue;
                
                k = single( (0:3) );
                tempB1 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*0 ) .* tempB(1:4), obj.EncodeMethod ) );
                tempB2 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*1 ) .* tempB(5:8), obj.EncodeMethod ) );
                tempB3 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*2 ) .* tempB(9:12), obj.EncodeMethod ) );
                tempB4 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*3 ) .* tempB(13:16), obj.EncodeMethod ) );
                tempB5 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*4 ) .* tempB(17:20), obj.EncodeMethod ) );
                tempB6 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*5 ) .* tempB(21:24), obj.EncodeMethod ) );
                tempB7 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*6 ) .* tempB(25:28), obj.EncodeMethod ) );
                tempB8 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*7 ) .* tempB(29:32), obj.EncodeMethod ) );
                tempB9 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*8 ) .* tempB(33:36), obj.EncodeMethod ) );
                tempB10 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*9 ) .* tempB(37:40), obj.EncodeMethod ) );
                tempB11 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*10 ) .* tempB(41:44), obj.EncodeMethod ) );
                tempB12 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*11 ) .* tempB(45:48), obj.EncodeMethod ) );
                tempB13 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*12 ) .* tempB(49:52), obj.EncodeMethod ) );
                tempB14 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*13 ) .* tempB(53:56), obj.EncodeMethod ) );
                tempB15 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*14 ) .* tempB(57:60), obj.EncodeMethod ) );
                tempB16 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*15 ) .* tempB(61:64), obj.EncodeMethod ) );

                tempB17 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*16 ) .* tempB(65:68), obj.EncodeMethod ) );
                tempB18 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*17 ) .* tempB(69:72), obj.EncodeMethod ) );
                tempB19 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*18 ) .* tempB(73:76), obj.EncodeMethod ) );
                tempB20 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*19 ) .* tempB(77:80), obj.EncodeMethod ) );
                tempB21 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*20 ) .* tempB(81:84), obj.EncodeMethod ) );
                tempB22 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*21 ) .* tempB(85:88), obj.EncodeMethod ) );
                tempB23 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*22 ) .* tempB(89:92), obj.EncodeMethod ) );
                tempB24 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*23 ) .* tempB(93:96), obj.EncodeMethod ) );
                tempB25 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*24 ) .* tempB(97:100), obj.EncodeMethod ) );
                tempB26 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*25 ) .* tempB(101:104), obj.EncodeMethod ) );
                tempB27 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*26 ) .* tempB(105:108), obj.EncodeMethod ) );
                tempB28 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*27 ) .* tempB(109:112), obj.EncodeMethod ) );
                tempB29 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*28 ) .* tempB(113:116), obj.EncodeMethod ) );
                tempB30 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*29 ) .* tempB(117:120), obj.EncodeMethod ) );
                tempB31 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*30 ) .* tempB(121:124), obj.EncodeMethod ) );
                tempB32 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*31 ) .* tempB(125:128), obj.EncodeMethod ) );
                
                tempB33 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*32 ) .* tempB(129:132), obj.EncodeMethod ) );
                tempB34 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*33 ) .* tempB(133:136), obj.EncodeMethod ) );
                tempB35 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*34 ) .* tempB(137:140), obj.EncodeMethod ) );
                tempB36 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*35 ) .* tempB(141:144), obj.EncodeMethod ) );
                tempB37 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*36 ) .* tempB(145:148), obj.EncodeMethod ) );
                tempB38 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*37 ) .* tempB(149:152), obj.EncodeMethod ) );
                tempB39 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*38 ) .* tempB(153:156), obj.EncodeMethod ) );
                tempB40 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*39 ) .* tempB(157:160), obj.EncodeMethod ) );
                tempB41 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*40 ) .* tempB(161:164), obj.EncodeMethod ) );
                tempB42 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*41 ) .* tempB(165:168), obj.EncodeMethod ) );
                tempB43 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*42 ) .* tempB(169:172), obj.EncodeMethod ) );
                tempB44 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*43 ) .* tempB(173:176), obj.EncodeMethod ) );
                tempB45 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*44 ) .* tempB(177:180), obj.EncodeMethod ) );
                tempB46 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*45 ) .* tempB(181:184), obj.EncodeMethod ) );
                tempB47 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*46 ) .* tempB(185:188), obj.EncodeMethod ) );
                tempB48 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*47 ) .* tempB(189:192), obj.EncodeMethod ) );
                
                tempB49 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*48 ) .* tempB(193:196), obj.EncodeMethod ) );
                tempB50 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*49 ) .* tempB(197:200), obj.EncodeMethod ) );
                tempB51 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*50 ) .* tempB(201:204), obj.EncodeMethod ) );
                tempB52 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*51 ) .* tempB(205:208), obj.EncodeMethod ) );
                tempB53 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*52 ) .* tempB(209:212), obj.EncodeMethod ) );
                tempB54 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*53 ) .* tempB(213:216), obj.EncodeMethod ) );
                tempB55 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*54 ) .* tempB(217:220), obj.EncodeMethod ) );
                tempB56 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*55 ) .* tempB(221:224), obj.EncodeMethod ) );
                tempB57 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*56 ) .* tempB(225:228), obj.EncodeMethod ) );
                tempB58 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*57 ) .* tempB(229:232), obj.EncodeMethod ) );
                tempB59 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*58 ) .* tempB(233:236), obj.EncodeMethod ) );
                tempB60 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*59 ) .* tempB(237:240), obj.EncodeMethod ) );
                tempB61 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*60 ) .* tempB(241:244), obj.EncodeMethod ) );
                tempB62 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*61 ) .* tempB(245:248), obj.EncodeMethod ) );
                tempB63 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*62 ) .* tempB(249:252), obj.EncodeMethod ) );
                tempB64 = fft( ComplexBlock32( exp( -1j*2*pi/N*k*63 ) .* tempB(253:256), obj.EncodeMethod ) );

                tempC = zeros( size( tempB ) );
                tempC(1:64:193) = tempB1.encodedValue;
                tempC(2:64:194) = tempB2.encodedValue;
                tempC(3:64:195) = tempB3.encodedValue;
                tempC(4:64:196) = tempB4.encodedValue;
                tempC(5:64:197) = tempB5.encodedValue;
                tempC(6:64:198) = tempB6.encodedValue;
                tempC(7:64:199) = tempB7.encodedValue;
                tempC(8:64:200) = tempB8.encodedValue;
                tempC(9:64:201) = tempB9.encodedValue;
                tempC(10:64:202) = tempB10.encodedValue;
                tempC(11:64:203) = tempB11.encodedValue;
                tempC(12:64:204) = tempB12.encodedValue;
                tempC(13:64:205) = tempB13.encodedValue;
                tempC(14:64:206) = tempB14.encodedValue;
                tempC(15:64:207) = tempB15.encodedValue;
                tempC(16:64:208) = tempB16.encodedValue;

                tempC(17:64:209) = tempB17.encodedValue;
                tempC(18:64:210) = tempB18.encodedValue;
                tempC(19:64:211) = tempB19.encodedValue;
                tempC(20:64:212) = tempB20.encodedValue;
                tempC(21:64:213) = tempB21.encodedValue;
                tempC(22:64:214) = tempB22.encodedValue;
                tempC(23:64:215) = tempB23.encodedValue;
                tempC(24:64:216) = tempB24.encodedValue;
                tempC(25:64:217) = tempB25.encodedValue;
                tempC(26:64:218) = tempB26.encodedValue;
                tempC(27:64:219) = tempB27.encodedValue;
                tempC(28:64:220) = tempB28.encodedValue;
                tempC(29:64:221) = tempB29.encodedValue;
                tempC(30:64:222) = tempB30.encodedValue;
                tempC(31:64:223) = tempB31.encodedValue;
                tempC(32:64:224) = tempB32.encodedValue;
                
                tempC(33:64:225) = tempB33.encodedValue;
                tempC(34:64:226) = tempB34.encodedValue;
                tempC(35:64:227) = tempB35.encodedValue;
                tempC(36:64:228) = tempB36.encodedValue;
                tempC(37:64:229) = tempB37.encodedValue;
                tempC(38:64:230) = tempB38.encodedValue;
                tempC(39:64:231) = tempB39.encodedValue;
                tempC(40:64:232) = tempB40.encodedValue;
                tempC(41:64:233) = tempB41.encodedValue;
                tempC(42:64:234) = tempB42.encodedValue;
                tempC(43:64:235) = tempB43.encodedValue;
                tempC(44:64:236) = tempB44.encodedValue;
                tempC(45:64:237) = tempB45.encodedValue;
                tempC(46:64:238) = tempB46.encodedValue;
                tempC(47:64:239) = tempB47.encodedValue;
                tempC(48:64:240) = tempB48.encodedValue;
                
                tempC(49:64:241) = tempB49.encodedValue;
                tempC(50:64:242) = tempB50.encodedValue;
                tempC(51:64:243) = tempB51.encodedValue;
                tempC(52:64:244) = tempB52.encodedValue;
                tempC(53:64:245) = tempB53.encodedValue;
                tempC(54:64:246) = tempB54.encodedValue;
                tempC(55:64:247) = tempB55.encodedValue;
                tempC(56:64:248) = tempB56.encodedValue;
                tempC(57:64:249) = tempB57.encodedValue;
                tempC(58:64:250) = tempB58.encodedValue;
                tempC(59:64:251) = tempB59.encodedValue;
                tempC(60:64:252) = tempB60.encodedValue;
                tempC(61:64:253) = tempB61.encodedValue;
                tempC(62:64:254) = tempB62.encodedValue;
                tempC(63:64:255) = tempB63.encodedValue;
                tempC(64:64:256) = tempB64.encodedValue;
                
                obj2 = ComplexBlock32( single( tempC ), obj.EncodeMethod );                
            else
                error('Error: input size not supported');
            end
        end
    end
end
