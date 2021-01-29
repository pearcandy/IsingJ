"""
  I s i n g J

  ver 0.2.1 Released Sep 2020 
  This code is distributed under the constitution of GNLv3.
  Log of ising.jl  
  2019/04/24  Released by PearCandy                 """

module IsingJ
    include("./Params.jl")
    include("./MCsampling.jl")
    include("./Initial.jl")
    include("./Smoothing_mag.jl")
    using ArgParse
    using Dates
    using Printf
    using Random
    using .Params
    using .MCsampling
    using .Initial
    using .Smoothing_mag
    export main

    # main
    function main()
        parser= ArgParseSettings("Julia code for solving 2D Ising model") # parser
        @add_arg_table parser begin
            "--lower_temp", "-l"
            arg_type = Float64
            default  = 0.1
            help     = "lowest temperature"        
            "--upper_temp", "-u"
            arg_type = Float64
            default  = 5.0
            help     = "highest temperature"        
            "--step_num", "-n"
            arg_type = Int
            default  = 30
            help     = "step number of temperature"        
            "--sweep_num", "-s"
            arg_type = Int
            default  = 2000
            help     = "montecarlo sweep number"        
            "--warm_up", "-w"
            arg_type = Int
            default  = 2000
            help     = "warm up steps before sampling"                
            "--exch_j", "-j"
            arg_type = Float64
            default  = 1.0
            help     = "exchange interaction"        
            "--ext_field", "-e"
            arg_type = Float64
            default  = 0.0
            help     = "magnetic field"        
            "--site_size", "-N"
            arg_type = Int
            default  = 64
            help     = "simulation size"
        end
            args = parse_args(parser)        
      
    """ Arguments        
        N_site    : number of site
        j0_val    : temperature
        low_t     : minimum temperature
        high_t    : maximum temperature
        ext       : external field
        inum_t    : number of temperature steps
        sweep_num : sweep number
        warm_up   : number of "warming-up" iteration  """

        N_site    = args["site_size"]
        j0_val    = args["exch_j"]
        low_t     = args["lower_temp"]
        up_t      = args["upper_temp"]
        ext_val   = args["ext_field"]
        inum_t    = args["step_num"]
        sweep_num = args["sweep_num"]
        warm_up   = args["warm_up"]
        cnst      = 2.26918  # Tc for square lattice 

        # set variables
        site_size=64
        exch_j=1.0
        lower_temp=0.1
        upper_temp=5.0
        ext_field=0.0
        step_num=30
        sweep_num=2000
        warm_up=2000

        params = Params.Variables(site_size,exch_j,lower_temp,upper_temp,
                                  ext_field,step_num,sweep_num,warm_up)
        N_site    = params.site_size
        j0_val    = params.exch_j
        low_t     = params.lower_temp
        up_t      = params.upper_temp
        ext_val   = params.ext_field
        inum_t    = params.step_num
        sweep_num = params.sweep_num
        warm_up   = params.warm_up

        # initialization of the physical quantities
        spin      = zeros(N_site,N_site)
        mag       = 0
        mag_ave   = 0
        dt        = (up_t-low_t)/float(inum_t)
        _mag      = []

        Random.seed!(100)
        spin_ini = cold_start(N_site)  # Initialization, random_start(N_site) 
        td=Dates.today()               # start motecalro simulation
        println(" ")
        println(td)
        println(" ")
        println("     I s i n g J")
        println("     Version 0.2.0, Released July 2019")
        println(" ")
        println("     Copyright (C) yassan1980")
        println(" ")
        println("     Usage: julia Ising.jl -h")
        println("---------------------------------------")        
        open("log.dat","w") do io
            println(io,"#temperature magnetization energy Cv Chi Tc")
            for i in 1:inum_t       # main routine (montecalro sampling) 
                _temp = low_t + dt*float(i)
                spin,mag,eng,cv,chi = mcsampling(spin_ini,sweep_num,warm_up,N_site,j0_val,ext_val,_temp)
                @printf(io, " %+06f"  ,_temp)
                @printf(io, " %+06f"  , mag)
                @printf(io, " %+06f"  , eng)
                @printf(io, " %+06f"  , cv)
                @printf(io, " %+06f"  , chi)
                @printf(io, " %+06f\n", cnst)
                append!(_mag,abs(mag))

                if i%10==0          # log 
                    println(i,"/",inum_t,"steps")
                else
                end
            end
        end
        
        output_T_vs_mag(low_t,dt,_mag)   # smoothing T vs magnetization

        #-- plot (using python) --
        #    pyplot=pyimport("matplotlib.pyplot")
        #    numpy=pyimport("numpy")
        #    x=[i for i in 1:N_site]
        #    y=[i for i in 1:N_site]
        #    xx,yy=numpy.meshgrid(x,y)
        #    pyplot.pcolormesh(xx,yy,spin,cmap="bwr")
        #    pyplot.axes().set_aspect("equal")
        #    pyplot.show()

    end 
end # module
#=-----------------------------------------------------------=#
#if occursin(PROGRAM_FILE,@__FILE__)
#    using .IsingJ
#    main()
#end

