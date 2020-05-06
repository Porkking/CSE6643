function [x] = nnls2(E,f)
% step 1:
[m,n]=size(E);
Z = 1:1:n;
P = [];
Universe = 1:1:n;
x = zeros(n,1);
% step 2
w = E'*(f-E*x);
while (true)
	% j in Z, w_Z
	w_Z = w(Z);
	% step 3
	if(isempty(Z) || max(w_Z) <=0)
		break;
	else
		% find t step 4
		t = Z(w_Z==max(w_Z));
		% step 5
		Z(Z==t) = [];
		P = setdiff(Universe, Z);
		% step 6
		while (true)
            E_p = zeros(m,length(P));
            for column= 1:length(P)
                E_p(:,column)= E(:,P(column));
            end
			% compute least square step 6
            %Solve full rank least square problem here!!!
			%z_test = E_p\f;
            %z_test=givens(E_p,f);
            %z_test=classicalgs(E_p,f);
            %z_test=householder(E_p,f);
            z_test=modifiedgs(E_p,f);
            %z_test=normaleqn(E_p,f);
            z = zeros(n,1);
            for i = 1:length(P)
                z(P(i)) = z_test(i);
            end
			% step 7
			z_P = z(P);
			if min(z_P) > 0
				break;
			else
				% step 8
				alphaArray = [];
				for element = P
					if z(element)<=0
						alphaArray = [alphaArray,x(element)/(x(element)-z(element))];
					end
				end
				% step 9
				alpha = min(alphaArray);
				% step 10
				x = x + alpha.*(z-x);
				for j = P
					if x(j) == 0
						Z = [Z,j];
					end
				end
				P = setdiff(Universe, Z);
			end
		end
		x = z;
		w = E'*(f-E*x);
	end
end
end