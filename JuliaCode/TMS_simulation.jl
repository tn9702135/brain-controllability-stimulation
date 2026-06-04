# Simulation of TMS in human cortical network

# Loading packages
using DifferentialEquations
using MAT
using DSP
using SparseArrays
using LinearAlgebra


# Setting parameters
const aee=1; const aie=1.4; const aei=0.7; const aii=0.4;
const kee=1.5; const kie=1; const kei=2; const kii=3;
const vee=10; const vie=8; const vei=-8; const vii=-12;
const de = 0.5; const di = 0.5;
const ni0=-20; const ne0=20;


const t_start_stim=100;
const t_during_stim=50;
const t_total=250;
const freq =1.0;
const peak =100;


# Load SC matrix
file = matopen("StructuralMatrix30.mat");
Structural_Matrix = read(file,"A1");
close(file)

# Load initial conditions
file = matopen("IntCon30.mat");
IntCon=read(file, "IntCon");
close(file);





N_sub=size(Structural_Matrix,2);
const Nn = size(Structural_Matrix[1],1) # Number of nodes

#println("Number of threads: $(Threads.nthreads())")
tid = Base.parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
for i =tid:tid
 C = 0.025*Structural_Matrix[i];
 IC = IntCon[i];
 # Adding self-connections
 
 C=C+kee*Diagonal((ones(Nn,Nn)))
 Nc=sum(!iszero,C); # Number of connections
 # Storing node degrees
 Cn_temp=zeros(size(C));
 Cn_temp[findall(!iszero,C)]=ones(size(findall(!iszero,C)));
 Cn=sum(Cn_temp,dims=1);
 # Storing number of ODEs for each node
 Vn=zeros(1,Nn+1);
  Vn[2:end]=cumsum(2*Cn+repeat([10],1,Nn),dims=2);
  # Integrating system
  tspan = (0.0,t_total);
  u0 = IC;
  Threads.@threads for j= 1:Nn
        target= j # Target region

        # Defining functions
#f(x,y)=(repeat([1],length(x),1)-x.^2-y.^2)./(repeat([1],length(x),1)+2*x+x.^2+y.^2)
f(x,y)=(1 .-x.^2-y.^2)./(1 .+2*x+x.^2+y.^2)
RF(x,y,n,d)=x.*y-y-d*(x.^2-y.^2+2*x+repeat([1],length(x),1))/2+n*(-x.*y-y);
IF(x,y,n,d)=(-x.^2+y.^2+2*x-repeat([1],length(x),1)-d*(2*x.*y+2*y)+n*(x.^2-y.^2+2*x+repeat([1],length(x),1)))/2;
RG(x,y,g,v)=sum((-x.*y-y)*v.*g-(x.^2-y.^2-repeat([1],length(x),1)).*g/2);
IG(x,y,g,v)=sum((x.^2-y.^2+2*x+repeat([1],length(x),1))*v.*g/2-x.*y.*g);

# TMS pulse function
function TMS_pulse(t,freq,peak)
    w=freq*1000; tau=0.00008;
    ts = mod(t,1/freq);
    if t>t_start_stim && t<(t_during_stim+t_start_stim)
        pulse = peak;
        #pulse=peak*sin(w*ts)*exp(-ts/tau);
        return pulse;
    else
        return 0;
    end
end




# Next-generation neural mass model
function CB_TMS(dydt,S,p,t)

    for n=1:Nn

        if n==target
            ne=ne0+TMS_pulse(t,freq,peak);
            ni=ni0+TMS_pulse(t,freq,peak);
        else
            ne=ne0;
            ni=ni0;
        end

        o=trunc(Int,Vn[n]);
        co=trunc(Int,Cn[n]);
        Fn=zeros(trunc(Int, Vn[n+1]-Vn[n]),1);
        xe=[S[o+1]]; ye=[S[o+2]]; xi=[S[o+3]]; yi=[S[o+4]];
        gii=[S[o+5]]; gie=[S[o+6]]; gei=[S[o+7]]; gee=S[o+8:o+7+co];
        hii=[S[o+8+co]]; hie=[S[o+9+co]]; hei=[S[o+10+co]]; hee=S[o+11+co:o+10+2*co];

        Fn[1]=RF(xe,ye,ne,de)[1]+RG(xe,ye,gei,vei)+RG(xe,ye,gee,vee);
        Fn[2]=IF(xe,ye,ne,de)[1]+IG(xe,ye,gei,vei)+IG(xe,ye,gee,vee);
        Fn[3]=RF(xi,yi,ni,di)[1]+RG(xi,yi,gii,vii)+RG(xi,yi,gie,vie);
        Fn[4]=IF(xi,yi,ni,di)[1]+IG(xi,yi,gii,vii)+IG(xi,yi,gie,vie);

        Fn[5:7+co]=[hii;hie;hei;hee];
        Fn[8+co]=dropdims(aii^2*(kii*f(xi,yi)-gii-2/aii*hii),dims=1)[1];
        Fn[9+co]=dropdims(aie^2*(kie*f(xe,ye)-gie-2/aie*hie),dims=1)[1];
        Fn[10+co]=dropdims(aei^2*(kei*f(xi,yi)-gei-2/aei*hei),dims=1)[1];
        Fn[11+co:10+2*co]=aee^2*(C[findall(!iszero,C[:,n]),n].*f(S[round.(Int,Vn[findall(!iszero,C[n,:])])+repeat([1],co,1)],S[round.(Int,Vn[findall(!iszero,C[n,:])])+repeat([2],co,1)])-gee-2/aee*hee);

        dydt[o+1:trunc(Int,Vn[n+1])]=Fn;
    end
end
         prob = ODEProblem(CB_TMS,u0,tspan);
         alg = MethodOfSteps(Tsit5());
          #alg = MethodOfSteps(RK4());
         #alg = Tsit5();
         alg = RK4() ;
         sol = solve(prob, maxiters = 1e8,alg, progress=true,save_idxs = append!(round.(Int,Vn[1:end-1].+1)[:,1],round.(Int,Vn[1:end-1].+2)[:,1]));


         # Save time series for firing rate of excitatory populations
         U=zeros(size(sol.u)[1],2*Nn);
         for n=1:size(sol.u)[1]
              U[n,:]=sol.u[n];
         end
          U=f(U[:,1:Nn],U[:,Nn+1:end]);

          # effect stimulation

         ind_t_before = findall(x->x<=t_start_stim-t_during_stim,sol.t);
         ind_t_start_stim = findall(x->x<=t_start_stim,sol.t);                  #time_start stimulation
         ind_t_stop_stim = findall(x->x<=t_start_stim+t_during_stim,sol.t);     #time_stop stimulation
         int_t_post_stim_delay_start = findall(x->x<=t_start_stim+2*t_during_stim,sol.t); # time_post_stimulation/start point
         int_t_post_stim_delay_end = findall(x->x<=t_start_stim+3*t_during_stim,sol.t); # time_post _stimulation/end point
         U_R_before_stim = U[ind_t_before[end]+1:ind_t_start_stim[end],:];      #firing rate before stimulation
         U_R_during_stim = U[ind_t_start_stim[end]+1:ind_t_stop_stim[end],:];   #firing rate during stimulation
         U_R_post_stim = U[int_t_post_stim_delay_start[end]+1:int_t_post_stim_delay_end[end],:]; #firing rate post stimulation


         # Compute FC before stimulation
         U_trans_before_stim =imag(hilbert(U_R_before_stim));
         R_before_stim=zeros(Nn,Nn);
          for n=1:Nn-1
              for m=n+1:Nn
                 R_before_stim[n,m]=abs((1/size(U_trans_before_stim)[1]*sum(exp.(im*(U_trans_before_stim[:,m]-U_trans_before_stim[:,n])))));
                R_before_stim[m,n]=R_before_stim[n,m];
               end
           end
 

          # Compute FC during stimulation
         U_trans_during_stim=imag(hilbert(U_R_during_stim));
         R_during_stim=zeros(Nn,Nn);
          for n=1:Nn-1
             for m=n+1:Nn
                 R_during_stim[n,m]=abs((1/size(U_trans_during_stim)[1]*sum(exp.(im*(U_trans_during_stim[:,m]-U_trans_during_stim[:,n])))));
                 R_during_stim[m,n]=R_during_stim[n,m];
                end
           end

         # Compute FC post stimulation
         U_trans_post_stim=imag(hilbert(U_R_post_stim));
         R_post_stim=zeros(Nn,Nn);
          for n=1:Nn-1
             for m=n+1:Nn
                 R_post_stim[n,m]=abs((1/size(U_trans_post_stim)[1]*sum(exp.(im*(U_trans_post_stim[:,m]-U_trans_post_stim[:,n])))));
                 R_post_stim[m,n]=R_post_stim[n,m];
               end
            end


     # Save time series and FC
     # file = matopen("timeseries_before_stim.mat", "w");
     # write(file, "Rf_before_stim", U_R_before_stim);
     # close(file);

      # file = matopen("timeseries_during_stim.mat", "w");
     # write(file, "Rf_during_stim", U_R_during_stim);
     # close(file);

     # file = matopen("timeseries_t.mat", "w");
     # write(file, "t", sol.t);
     # close(file);

     file = matopen("FC_before_stim$i$j.mat", "w");
     write(file, "R_before_stim", R_before_stim);
     close(file);

     file = matopen("FC_during_stim$i$j.mat", "w");
     write(file, "R_during_stim",R_during_stim);
     close(file);

     file = matopen("FC_post_stim$i$j.mat", "w");
     write(file, "R_post_stim",R_post_stim);
     close(file);
  end
  

# Save solution and time.
# file = matopen("steadystate.mat", "w");
# write(file, "sol", sol.u);
# close(file);

# file = matopen("steadystate_time.mat", "w");
# write(file, "t", sol.t);
# close(file);
end
