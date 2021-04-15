function expr = Model(t,in2,in3)
%MODEL
%    EXPR = MODEL(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    14-Apr-2021 09:36:20

He = in2(6,:);
Hi = in2(2,:);
Hr = in2(4,:);
Hs = in2(5,:);
Me = in2(3,:);
Mi = in2(8,:);
Ms = in2(7,:);
param1 = in3(:,1);
param2 = in3(:,2);
param3 = in3(:,3);
param4 = in3(:,4);
param5 = in3(:,5);
param6 = in3(:,6);
param7 = in3(:,7);
t2 = He.*param6;
t3 = Hi.*param3;
t4 = Me+Mi+Ms;
t5 = Hi.*Ms.*param2;
t6 = He.*Me.*2.0;
t7 = He.*Mi.*2.0;
t8 = He.*Ms.*2.0;
t9 = Hs.*Mi.*param1.*5.0;
t10 = He+Hi+Hr+Hs;
t11 = -t9;
t12 = 1.0./t4;
t13 = 1.0./t10;
expr = [t2;Hi.*(-2.0./5.0)+t2-t3;-t13.*(-t5+He.*Me.*param5+He.*Me.*param7+Hi.*Me.*param5+Hi.*Me.*param7+Hr.*Me.*param5+Hr.*Me.*param7+Hs.*Me.*param5+Hs.*Me.*param7);Hr.*(-2.0./5.0)+t3;(t12.*(t6+t7+t8+t11+Hi.*Me.*2.0+Hi.*Mi.*2.0+Hr.*Me.*2.0+Hr.*Mi.*2.0+Hi.*Ms.*2.0+Hr.*Ms.*2.0))./5.0;t12.*(t6+t7+t8+t11+Me.*t2.*5.0+Mi.*t2.*5.0+Ms.*t2.*5.0).*(-1.0./5.0);-t13.*(t5-He.*param4-Hi.*param4-Hr.*param4-Hs.*param4+He.*Ms.*param5+Hi.*Ms.*param5+Hr.*Ms.*param5+Hs.*Ms.*param5);Me.*param7-Mi.*param5];