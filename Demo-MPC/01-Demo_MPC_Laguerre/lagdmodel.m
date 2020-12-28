clear
numd=[1 -0.1];
dend=conv([1 -0.8],[1 -0.9]);
N_sim=60;
k=1:(N_sim);
H=dimpulse(numd,dend,k);

a=0;
N=4;
[A1,L0]=lagd(a,N);
L(:,1)=L0;
for kk=2:N_sim
L(:,kk)=A1*L(:,kk-1);
end

c1=L(1,:)*H;
c2=L(2,:)*H;
c3=L(3,:)*H;
c4=L(4,:)*H;
H_model=c1*L(1,:)+c2*L(2,:)+c3*L(3,:)+c4*L(4,:);
figure (2)
plot(k,H)
hold on
plot(k,H_model,'LineWidth',2,'Color',[.8 0 0]);
set(gca,'FontSize',20,'FontName','helvetica');
legend('data','model')
xlabel('Sampling Instant');
ylabel('Impulse Response')
