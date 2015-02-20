%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Matlab file             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
%%%%            Last modified on 07/26/2005          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

hh = 2;  %profondita di partenza%
dh = 0;
t  = 0;
dt = 1800; %1sec

C=0; %concentrazione di partenza

L=0;    
rand('state',sum(100*clock));
for j = 1 : 2000,

   j;
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% DURATA VENTO  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   f = rand;     % f is a pseudo random number

   x = ((1.0975-f)/42.64)^(-1.86);
%   x = 687.23 * (f^(-0.6504));   %secondi
   durata = x - mod(x,dt) + dt;
   

   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% INTENSITA' VENTO  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   d = abs(rand);   
   
   % if (-2.4258 * log(d*1000) + 17.57) < 40    
     % uw = -2.4258 * log(d*1000) + 17.57;
     % else
     % uw = 40;
     % end
     
     uw = -5.62 * log((1.0118-d)/1.579);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   CICLO MAREA     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   for jj = t : dt : durata+t,

%      eta = 0.35 * sin(2*pi*(jj)/43200);

% tidal elevation 


   eta=0.254*sin(2*pi/360*(28.98*jj/3600+120))+0.143*sin(2*pi/360*(30*jj/3600+324))+0.042*sin(2*pi/360*(28.43*jj/3600+343))+0.032*sin(2*pi/360*(30.08*jj/3600+117))+0.155*sin(2*pi/360*(15.04*jj/3600+82))+0.04*sin(2*pi/360*(13.94*jj/3600+256))+0.057*sin(2*pi/360*(14.96*jj/3600+97))+0.015*sin(2*pi/360*(15*jj/3600+262));
  

h = hh + eta;

      if h < 0
         
          
          t = t + dt;
            L=L+1;
          %%%%%%%%%% se si supera il livello del mare,non
                           %%%%%%%%%% si ha piu formazione di
                           %%%%%%%%%% onde------->bisogna pero aggiornare i
                           %%%%%%%%%% dati
                      
                            
         
      vento(L) = uw;
      altezza(L) = hh;
      tempo(L) = t;
      corrente(L) = 0;
      livello(L) = eta;
      onda(L) = 0;
      taumo(L) = 0;
      tauidro(L) = 0;
      tautot(L) = 0;
      
      
 %     PLOTTA
      
           break         
           
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  MODELLO ONDE     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
       
      c  = 0.015;       %coeff per bottom friction%
      cd = 0.0012;
      Hmax = 0.78 * h;  %per avere breaking%
      g = 9.806;
      Rhow = 1027;      %densita H2o%
      r = 0.001217137;  %dens airia/densita H2o%
      T = 2;            %periodo onda%
      apm = 0.00457;    %coeff per whitecapping%
      Dt = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  NUMERO D'ONDA    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      kk = 0.5;
      scarto = 1;
      T = 2;
      w = 2*pi/T;

      cit=0;

      while ((scarto>0.001) & (cit<100)),
         cit = cit + 1;
         fk = sqrt(g*kk*tanh(kk*h)) - w;
         dfk=0.5*sqrt(g/(kk*tanh(kk*h))*(tanh(kk*h)-kk*h*tanh(kk*h)^2+kk*h));
         Df = -fk/dfk;
         k = kk + Df;
         scarto = abs(k-kk);
         kk=k;
      end

      %---------------------->k deve cambiare con le h (depht)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  MARE COMPLETAMENTE SVILUPPATO      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

      N=0;
      toll=1;
      H=0;
      nit=0;

      while ((toll>0.00000029) & (nit<500));

         nit=nit+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  white capping                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         awc = (1/8*H^2) * w^(4/g^2);
         Wwc = -0.000033 * w * (awc/apm)^2;
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%    bottom friction                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         Wbf = (-4*c*pi*H*k)/(T*sinh(k*h)*sinh(2*k*h));
         
         b = 5*r/T*(uw*k/w-0.9);
         a = ((80*r^2*w*cd^2*uw^4)/((g*k)^2))/w;
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%           breaking                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
         
         
         Qb = 0;

         if H/Hmax < 0.4
            Qb = 0;
         elseif H/Hmax<1
            Qb = 1.8*(H/Hmax-0.4)^3+1.7*(H/Hmax-0.4)^2;
         else
            Qb=1;
         end

         Wbr = -2*Qb/T*(H/Hmax)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%     bilancio energia dissipata      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

         NN=N+(a+b*N+Wbf*N+Wwc*N+Wbr*N)*Dt;
         if(NN<0)
             NN=0;
         end
         HH=sqrt((8*(NN*w)/(Rhow*g)));
         toll = abs(HH-H);
         H = HH;
         N = NN;
      end

      Hfas = HH;

      Uwf = pi * Hfas/(T*sinh(k*h));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%    dati                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

      D50 = 0.00002;
      rhoD50 = 2650;
      ss = 2.580331061;
      visc = 0.00000101;
      
      r=1;
      vs=0.00032;
      
   %     qw=0.15;
      
    %  if hh>0
     %   UUU=qw/hh;
     % else
   %   UUU=0
   %end
      
 % if UUU>0.25
  %    UUU=0.25;              
           
  %end
      
%      UUU=0.05+0.045*hh;
    %UUU=0.15;
    
    UUU=-0.0086*hh^3+0.0198*hh^2+0.0513*hh+0.0625;
    
   %   Uc=abs(UUU*cos(2*pi*(jj)/43200));
      
 Uc=UUU*(0.344*cos(2*pi/360*(28.98*jj/3600+120))+0.194*cos(2*pi/360*(30*jj/3600+324))+0.057*cos(2*pi/360*(28.43*jj/3600+343))+0.043*cos(2*pi/360*(30.08*jj/3600+117))+0.21*cos(2*pi/360*(15.04*jj/3600+82))+0.054*cos(2*pi/360*(13.94*jj/3600+256))+0.077*cos(2*pi/360*(14.96*jj/3600+97))+0.02*cos(2*pi/360*(15*jj/3600+262)));
      
      Ks = 30;

     tidro= Rhow*g*Uc^2/(Ks^2);  % tidro= Rhow*g*qw^2/(Ks^2*h^(7/3)); % CARMEN, CHECK This line.


      Z = D50/12;
      A = Uwf*T/(2*pi);
      Re = Uwf*A/visc;
      fwr = 1.39*(A/Z)^(-0.52);

      if (Re <= 500000)
        fws = 2*Re^(-0.5);
      else
        fws = 0.0521*Re^(-0.187);
      end

      fws = fws;
      f   = max(fwr,fws);
      tmo = 0.5*Rhow*Uwf^2*f;

      ttot = tmo + tidro*(1+1.2*(tmo/(tidro+tmo))^(3.2));

      %if mod(floor(t/15552000),2)==0
      
      
       %       tcr = 1.1;
        %      X = 0.0029;
         %     b = 0.99;
         % else
              
           %   tcr = 0.7;
            %  X = 0.0012;
             % b = 2.13;
              %  end
      tcr=0.7;
      
      
      % erosion
      
      if (ttot > tcr)
    
%     Em = X *ttot^(b);      %kg/m2/sec%
      Em=4.12*0.0001*(ttot-tcr);
      
      else
         Em = 0;
      end

      % deposition
      
      if (ttot < 0.65)
      %  D = 0.87*10^(-4); %*/(h+0.35)/1.35;% 7*10^(-4);               %kg/m2/sec%
      D=r*vs*C;
      else
        D=0;
      end

      % change in elevation
      dh =(D-Em)*dt/(1800)*50;     %rho sed 1800 e porosita 0,1%
      
      % change in concentration 
      C=C+1/hh*(Em-D)*dt;
 
      
      % concentration of sediments in the river kg/m3
      Criver=0.02;
      % spatial dimension of the cell
      dx=100;
      
      % add sediment during ebb from a river 
      if Uc>0
          C=C+Uc*Criver*dt/dx;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
      
      % eliminate sediment during ebb
          C=C-0.2*Uc*C*dt/dx;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
      
      end
      
      hh = hh - dh;
      if mod(t,86400)==0
        hh
        j 
    end
    
      L = L + 1;
      t = t + dt;

      % save the variables in a matrix
      vento(L) = uw;
      altezza(L) = hh;
      tempo(L) = t;
      corrente(L) = Uc;
      livello(L) = eta;
      onda(L) = Hfas;
      taumo(L) = tmo;
      tauidro(L) = tidro;
      tautot(L) = ttot;
      concentrazione(L)=C;       
      if Uc>0
        input(L)=Uc*Criver*dt/dx;
        output(L)=0.2*Uc*C*dt/dx;
    else
        input(L)=0;
        output(L)=0;
    end
        
   end   % End of the for jj = 1 : dt : durata Loop.

end
%plot (tempo, altezza, 'k-');   % Plot stuff.

PLOTTA

save whosANNO

fidanno = fopen('anno','w'); % opens a file called data1.%
for kkk = 1:1:L,
  fprintf(fidanno,'%18.7f\n',altezza(kkk));  %print as a column.%
end,

fclose(fidanno);  % close the file with fid1.%






