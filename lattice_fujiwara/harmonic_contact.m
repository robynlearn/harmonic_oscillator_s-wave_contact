function [output] = harmonic_contact
% Energy is measured in hbar*omega
%
% scattering length is measured in terms of the reduced mass harmonic
% oscillator length a_ho = sqrt(hbar/(mu*omega))

% Energy to inverse scattering length
E2ainv = @(E) (2*gamma(-E/2+3/4)./gamma(-E/2+1/4));

% dE/d(1/a)
dEd1a = @(E) (gamma(1/4-E/2)./gamma(3/4-E/2))./(psi_custom(1/4-E/2)-psi_custom(3/4-E/2));

% Energy bounds for each branch
Ebounds = [-100 0.5:2:10.5]+1;
    
Np=1e5;

    output=struct;
    for jj=1:length(Ebounds)-1
        Evec=linspace(Ebounds(jj)+1E-6,Ebounds(jj+1)-1E-6,Np);
        Cvec = dEd1a(Evec);
        a_inv = E2ainv(Evec);       
        output(jj).E2dEd1a = @(Ein) interp1(Evec,Cvec,Ein);
        output(jj).ainv2dEd1a = @(a_inv_in) interp1(a_inv,Cvec,a_inv_in);    
        output(jj).ainv2E = @(a_inv_in) interp1(a_inv,Evec,a_inv_in);    
    end

end
