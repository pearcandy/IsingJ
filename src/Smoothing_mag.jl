module Smoothing_mag

    using Printf
    export output_T_vs_mag

    """
      plot_T_vs_mag(_temp,_mag)
      
      output T vs Magnetization
      _temp: temperature
      _mag : absolute values of magnetization

    """
    function output_T_vs_mag(low_t,dt,_mag)
        open("t_mag.dat","w") do io
            println(io,"# temperature vs ave.magnetization")

            for i in 1:length(_mag)
                _temp = low_t + dt*float(i)

                # smoothing of magnetization
                if(3<i<length(_mag)-3)
                    mag_ave=(_mag[i-2]+_mag[i-1]+_mag[i]+_mag[i+1]+_mag[i+2])/5
                else
                    mag_ave=_mag[i]
                end

                # -- output --
                @printf(io, " %+06f",     _temp  )
                @printf(io, " %+06f\n",   mag_ave)
            end
        end
    end # function

end #module
