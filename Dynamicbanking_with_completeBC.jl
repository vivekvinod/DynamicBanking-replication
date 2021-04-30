using Plots
using EconPDEs

struct BoltonLiWangYangModel
  r::Float64
  ρ::Float64
  αA::Float64
  σA::Float64
  δX::Float64
  σX::Float64
  θ::Float64
  ω::Float64
  ϕ::Float64
  ψ0 ::Float64
  ψ1::Float64  
  ξL::Float64
  ξK::Float64  
end


function BoltonLiWangYangModel(r = 0.01, ρ = 0.045, αA = 0.002, σA = 0.1, δX = 0.0, σX = 0.05, θ = 0.5, ω = 5.3, ϕ = 0.8, ψ0 = 0.001  , ψ1 = 0.05 ,ξL = 20.0,  ξK = 14.3)
 BoltonLiWangYangModel(r, ρ, αA, σA, δX, σX, θ, ω, ϕ, ψ0, ψ1, ξL, ξK)
end



klower = 1/(1-0.05)-1

# Only getting output for n = 2

# What are the bounds for k? the min bound is decided by the SLR requirement or capital requirement

function initialize_stategrid(m::BoltonLiWangYangModel; n = 400)
  OrderedDict(:k => range(klower,stop = 0.3, length = n))
end


function initialize_y(m::BoltonLiWangYangModel, stategrid::OrderedDict)
   OrderedDict(:j => stategrid[:k])
end



function (m::BoltonLiWangYangModel)(state::NamedTuple, y::NamedTuple)
  r = m.r; ρ = m.ρ; αA = m.αA;  σA = m.σA; δX = m.δX; σX = m.σX; θ = m.θ; ω = m.ω; ϕ = m.ϕ; ψ0 = m.ψ0; ψ1 = m.ψ1;ξL = m.ξL; ξK=m.ξK 
  k = state.k
  j, jk, jkk = y.j, y.jk, y.jkk
    
    
  #Taking log transformation of the state variable to avoid function -> inf because of second derivative in denominatore 
    # v=j
    #  vk=jk
    # vkk=jkk
    
  v = exp.(j)
  vk = exp.(j) .* jk
   vkk = exp.(j) .* (jk .* jk + jkk)
    

#Defining risk tolerance gamma    
  gamma_k = -vkk*k/vk
    
    
 #    Boundary of gamma_k which is obtained after solving the capital ratio of 14.3
    if  (gamma_k <  2/139)
      vkk = vk * (-2/139)/k 
   end
    
     gamma_k = -vkk*k/vk
            
      #Interest  
        
      i = (((v - vk * k )/ vk) - (1/ω))/(θ*ω)


 #  gamma_k = max(gamma_k, 2/139)
    
 pia = (k/(1+k))*(αA/(gamma_k *σA*σA ) + (σX * ϕ / σA))  
    
    
    
   vt =  (v - vk * k)*(ω*i - δX ) + (0.5 * vkk * k^2 * σX^2 ) + (vk)*(1+k)*(r+ pia*αA ) + 0.5*vkk*(1+k)^2 * ( pia*αA )^2 - vk*(i+0.5*θ*(ω*i)^2) - vkk*k*(1+k)*pia*σA*σX*ϕ - ρ*v 


    
  μw =   -(ω*i - δX )*k +(1+k)*(r+ pia*αA)-(i+0.5*θ*(ω*i)^2)


    
  return (vt,), (μw,), (j = j, jk = jk, jkk = jkk, k= k, pia = pia)
end



m = BoltonLiWangYangModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
y, result, distance = pdesolve(m, stategrid, y0; bc = OrderedDict(:jk=> (1.5, 0.5)))



#========================================================================================

Iterate on boundary conditions until the solution satisfies the conditions given in Bolton Et. al (2020)

========================================================================================#

function f(m, x, stategrid, y)
  y, result, distance = pdesolve(m, stategrid, y; bc = OrderedDict(:jk => (x[1], x[2])))
  u1, u2, k = result[:j], result[:jk], result[:k]
  out = zeros(4)      
     v = exp.(u1)
    vk = exp.(u1).*u2
    
    
  mi = argmin(abs.(1 + m.ψ1 .- vk))
    
    
    #Eqn 37 , getting the optimal issuance m (mi)
  out[1] = vk[mi] - 1 - m.ψ1
    
    #Eqn36
    
  out[2] =  v[mi] - m.ψ0 - (1 + m.ψ1) * k[mi] - v[1]
    
  wi = argmin(abs.(1 .-vk))
    
    #Eqn34 (Setting the FOC at the upper bound as 1 , )
  out[3] = vk[wi] - 1
    
    i = (v[wi] - k[wi] - 1/(m.ω) ) / (m.θ * m.ω)    
    
    gamma_k = -k[wi]
    
     pia = (k[wi]/(1+k[wi]))*(m.αA/(gamma_k *m.σA*m.σA ) + (m.σX * m.ϕ / m.σA))  
    
    #BC at the u[]
    
   out[4] =  (v[wi] - v[wi] * k[wi])*(m.ω*i - m.δX ) + (0.5 * v[wi] * k[wi]^2 * m.σX^2 ) + (v[wi])*(1+k[wi])*(m.r+ pia*m.αA ) + 0.5*v[wi]*(1+k[wi])^2 * ( pia*m.αA )^2 - v[wi]*(i+0.5*m.θ*(m.ω*i)^2) - v[wi]*k[wi]*(1+k[wi])*pia*m.σA*m.σX*m.ϕ - m.ρ * v[wi] 

    
 println(wi)
    
  #To add eqn  (Value of the ODE at the bound)
    
  # out[4] = m.ρ * u1[wi] -   (-(m.ω*i - m.δX )*k +(1+k)*(r+ pia*αA)-(i+0.5*θ*(ω*i)^2) )

    
  return out
end



using LeastSquaresOptim
newsol = optimize(x -> f(m, x, stategrid, y0), [1.0,  1.0], Dogleg())
# Solve the ODE using boundary conditions obtained from the optimizer

y, result, distance = pdesolve(m, stategrid, y0; bc = OrderedDict(:jk => (newsol.minimizer[1],newsol.minimizer[2])))

    


