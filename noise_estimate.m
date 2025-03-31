function [N1hat,N2hat]=noise_estimate(u,delta,L)
% FUNCTION [N1hat,N2hat]=noise1_estimate(u,delta,L)
%
% Produces vectors of noise estimates for first and second derivatives
% u: input signal
% delta: time step
% L: bound on the third derivative
N1hat=0;
N2hat=0;
    for ii=3:3:length(u)-1 % Index T
        for sig=3:3:ii     % Index sigma
            Q1=4*u(end)-(9/2)*u(end-(sig/3))+(1/2)*u(end-sig) - ...
                (4*u(end)-(9/2)*u(end-(ii/3))+(1/2)*u(end-ii))*sig/ii;
            N1upd=(abs(Q1)-L*sig*(ii^2-sig^2)*delta^3/18)/(9+sig/ii);
            N1hat=max(N1hat,N1upd);

            Q2=6*u(end)-9*u(end-(sig/3))+3*u(end-sig) - ...
                (6*u(end)-9*u(end-(ii/3))+3*u(end-ii))*(sig/ii)^2;
            N2upd=(abs(Q2)-4*L*sig^2*(ii-sig)*delta^3/9)/(18+6*(sig/ii)^2);
            N2hat=max(N2hat,N2upd);
        end
    end
