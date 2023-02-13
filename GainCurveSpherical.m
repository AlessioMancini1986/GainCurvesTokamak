%% Tokamak Gain Curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Point of contact: Alessio Mancini

%      e-mail= amancini@us.es

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%Parameter 

c=0.0491; %scaling factor

H=1 %H-factor

A=3.22; %Aspect ratio conventional tokamak SPARC
%A=2; %Aspect ratio spherical tokamaka
%A=3.22; SPARC Parameter
k=1.97; %Elongation conventional tokamak
%k=2.2; %Elongation spherical tokamak
%k=1.97; SPARC Parameter
d=0.5; %Triangularity
q=3.05; %here q is the q* SPARC Parameter

%% Plot of Q as a function of R,BT and IP as a function of R,BT for a q95>3

map=colormap2();
%Ip(k,d,A,map);
%hold on;
Pfus(q,A,map);
hold on;
Q_gain(c,H,k,q,A);
hold on

% Function for the plot of the plasma current
function []=Ip(k,d,A,map)
    % here the q95 is different from the q defined above which is the q*.
    % The important is that is greater than 2-3.
    %q=4.5;
    q=5.5;
    K=sqrt(0.5*(1+k));
    fd=(1+(k^2)*(1+2*d^2-1.2*d^3))/(1+k^2);
    fA=(1.17-0.65/A)/(1-1/(A^2))^2;
    HH=(5*K*K*fA*fd)/A/A/q;
    R=linspace(0,10,1000);
    BT=linspace(0,14,1000);

    for kk=1:length(BT)
        for ww=1:length(R)
        
            Ip(kk,ww)=HH*BT(kk)*R(ww);
            %RR(kk,ww)=R(ww);
            %BTT(kk,ww)=BT(kk);

        end
    end
    [BTT,RR]=meshgrid(BT,R);
    contourf(BTT,RR,Ip,100,'LineStyle','none');
    colorbar
    xlim([0, 14]); ylim([0 10]);
    caxis([0 20])
    hcb=colorbar;
    hcb.Title.String = "I_{P} [MA]";
    set(hcb,'Ticks',[0:2:30])
    colormap(map);
    hold on;
end

% Function for the plot of the fusion power 
function []=Pfus(q,A,map)
    % Scaling of the parameter to ITER tokamak

        P=400*1e+6; %MW
        q1=3.0;
        R1=6.2;
        A1=6.2/2.0;
        B1=5.3;
        betaN=1.05;
        cfus=P*(q1)^2*(A1)^4/(betaN)^2/(B1)^4/(R1)^3;
        R=linspace(0,10,100);
        BT=linspace(0,14,100);
        %A=1/0.31;
        bN=betaN;
        %q=3.05;
        cf=cfus;

            for i=1:length(BT)
                for j=1:length(R)
                   
%                     RR(i,j)=R(j);
%                     BTT(i,j)=BT(i);
                    PFUS(i,j)=cf*(bN)^2*((BT(i))^4*(R(j))^3)/(q)^2/(A)^4/1e+6;
                    
                end
            end

        [BTT,RR]=meshgrid(BT,R);
        %levels=0:1000/1000:10000/1000;
        levels=[0:1000:10000]/1000;
        %contourf(BTT,RR,PFUS/1000,levels,"LineColor","none");
        mesh(BTT,RR,PFUS/1000);
        xlabel('BT');ylabel('R');
        colorbar;
        %xlim([0, 14]); ylim([0 20]);
        %caxis([0 10])
        %hcb=colorbar;
        %hcb.Title.String = "P_{fus} [GW]";
        %set(hcb,'Ticks',levels)
        colormap(map);
        hold on;



end
% Function for the evaluation of the Q gain
function []=Q_gain(c,H,k,q,A)

r=linspace(0.01,15,1000);
Q=[0.1,1,2,10,50];
fn=strings(1,length(Q));
fn2=strings(1,length(c));

for j=1:length(c)     
    
    fn2(1,j)=strcat('c=',num2str(c(j)));
    D=5*c(j)*(H^2)*(k^3.5);
    E=5*(q^3)*A^3;
    F=(c(j)*(H^2)*(k^(3.5)));
    
    for i=1:length(Q)
    
            fn(1,i)=strcat('Q=',num2str(Q(i))); 
    
            S=((Q(i)*E)/(Q(i)*F + D)); 
            Bt=(S./(r.^2)).^(1/3);
            figure(j)
            plot(Bt,r,'LineWidth',2);
            xlabel('B_T [T]','FontSize',22); ylabel('R [m]','FontSize',22);
            xlim([0, 14]); ylim([0 10]);
            hold on;
            

    end

    % Plot of the different machine. It can be commented depending whether
    % you want to see the different position of the tokamaks on the plot.
%eSMART
%     plot(5,1.9,'bp','LineWidth',4);
%     text(4,6,'eSMART, I_{P}=8.2 MA, B_{T}=8 T, R_{0}=0.9 m, T_{0}=16 keV','Fontsize',14,'Color', 'k');
%     annotation('textbox','String', "Q=0.1",'LineStyle',"none",'FontSize',12,'Color','w');
%     annotation('textbox','String', "Q=1",'LineStyle',"none",'FontSize',12);
%     annotation('textbox','String', "Q=2",'LineStyle',"none",'FontSize',12);
%     annotation('textbox','String', "Q=10",'LineStyle',"none",'FontSize',12);
%     annotation('textbox','String', "Q=50",'LineStyle',"none",'FontSize',12);
%     hold on;
     
%     %SPARC
%     plot(12.2,1.85,'rp','LineWidth',7);
%     text(12.2,1.85,'Sparc','Fontsize',15,'Color', 'w');
%     hold on;
%     %CMOD
%     plot(8,0.67,'gx','LineWidth',7);
%     text(8,0.67,'C-mod','Fontsize',15,'Color', 'w');
%     hold on;
%     %AUG
%     plot(3.9,1.65,'gx','LineWidth',7);
%     text(3.9,1.65,'AUG','Fontsize',15,'Color', 'w');
%     hold on;
%     %DIIID
%     plot(2.2,1.66,'gx','LineWidth',7);
%     text(2.2,1.66,'DIIID','Fontsize',15,'Color', 'w');
%     hold on;
%     %JET
%     plot(3.45,2.96,'gx','LineWidth',7);
%     text(3.45,2.96,'JET','Fontsize',15,'Color', 'w');
%     hold on;
%     %IGNITOR
%     plot(13,1.32,'yd','LineWidth',7);
%     text(13,1.32,'Ignitor','Fontsize',15,'Color', 'w');
%     hold on;
%     %ARC
%     plot(9.2,3.3,'rp','LineWidth',7);
%     text(9.2,3.3,'ARC','Fontsize',15,'Color', 'w');
%     hold on;
% 
%     %cit
%     plot(10,2.1,'yd','LineWidth',7);
%     text(10,2.1,'CIT','Fontsize',15,'Color', 'w');
%     hold on;
%     %FIRE
%     plot(10,2.14,'yd','LineWidth',7);
%     text(10,2.14,'FIRE','Fontsize',15,'Color', 'w');
%     hold on;
%     %BPX
%     plot(9,2.59,'yd','LineWidth',7);
%     text(9,2.59,'BPX','Fontsize',15,'Color', 'w');
%     hold on;
%     %ITER
%     plot(5.3,6.2,'bs','LineWidth',7);
%     text(5.3,6.2,'ITER','Fontsize',15,'Color', 'w');
%     hold on;
%     %EU-DEMO
%     plot(5.9,9,'bs','LineWidth',7);
%     text(5.9,9,'EU-DEMO','Fontsize',15,'Color', 'w');
%     hold on;
    %lgd=legend(fn,'Fontsize',16);
    %legend show
    %title(fn2(1,j),'FontSize',16);
    %title('P_{fus}=200 MW','FontSize',16)
    set(gca,'FontSize',16);
%     hh(j)=gcf;
%     file_name=strcat(fn2(1,j),'_A=2.2_k=2.2.jpeg');
%     %nameimage=strcat(Folder,file_name);
%     saveas(hh(j),file_name);
%     %close;    
end

end
% Function for the PSFT Colormap 
function [map]=colormap2()

%%Custom color map: IDL Omega-II colormap

T = [0,   0,   0          %// dark
     0, 0,  255         %// blue
     255, 0, 0        %// red
     255, 255, 0        %// yellow
     255, 255, 255]./255; %// white 
 %Setting color intervals length
 x = [0
     70
     130
     200
     255];
 %Interpolation between colors
 map = interp1(x/255,T,linspace(0,1,255));
 %Color bar limits form 0 to 0.7 (black to white)
%  I = linspace(0,1,255);
% imagesc(I(ones(1,10),:)')
% colormap(map)
end
