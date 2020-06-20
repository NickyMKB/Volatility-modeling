function [NIC] = NICAS(alpha_pos,alpha_neg,residual)
T= size(residual,1);
% Extract the sample size (make sure returns are a column vector)
% Run the dynamic scale
for t=1:(T)
    if t==1;
        NIC(t,1)=0;
    else
        if (residual(t,1))>=0;
            %inidicator functions returns 1
            indicator =1;
        else
            indicator = 0;
        end
    NIC(t,1) = alpha_pos*(indicator*residual(t,1)^2 - 1/2) + alpha_neg*((1-indicator)*residual(t,1)^2 - 1/2);
end

end
% Close the function
end