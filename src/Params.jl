module Params
    mutable struct Variables
        site_size::Int64
        exch_j::Float64
        lower_temp::Float64
        upper_temp::Float64
        ext_field::Float64
        step_num::Int64
        sweep_num::Int64
        warm_up::Int64
    end
end

 """ Arguments        
     N_site    : number of site
     j0_val    : temperature
     low_t     : minimum temperature
     high_t    : maximum temperature
     ext       : external field
     inum_t    : number of temperature steps
     sweep_num : sweep number
     warm_up   : number of "warming-up" iteration  """

site_size=64
exch_j=1.0
lower_temp=0.1
upper_temp=5.0
ext_field=0.0
step_num=30
sweep_num=2000
warm_up=2000

#using .Params
#params=Params.Variables(site_size,exch_j,lower_temp,upper_temp,
#                        ext_field,step_num,sweep_num,warm_up)

