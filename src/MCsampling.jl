module MCsampling

    export mcsampling
    """
    Algorithms for solving classical Ising model   
    Args:                                                       
        spin:
            spin    
        sweepnum:
            toal  steps for sampling
        warm_up:
            warm up steps before sampling
        N:
            simulation (Matrix) size             
        j0:
            exchange coupling
        ext:
            magnetic field
        temp:
            temperature (Tc=2.26918)
    """
    function mcsampling(spin,sweepnum,warmup,N,j0,ext,temp)
        # physical quantities
        mag=0
        eng=0
        mag2=0
        eng2=0
        # montecarlo sweep
        for io in 1:sweepnum+warmup
            for i in 1:N
                for j in 1:N
                    E_old=get_energy(spin,i,j,j0,N,ext)
                    # spin flip at j-site
                    spin[i,j]=(-1.0)*spin[i,j]
                    E_new=get_energy(spin,i,j,j0,N,ext)
                    # energy diff after spin flip
                    delta=E_new-E_old
                    # judgement of spin flip or not
                    if(delta <= 0)
                        spin[i,j]=spin[i,j]
                    else
                        if(rand() <= exp((-1.0)*delta/temp))
                            spin[i,j]=spin[i,j]
                        else
                            spin[i,j]=(-1.0)*spin[i,j]
                        end                        
                    end
                end
            end

            if(io > warmup)
                # calculate phyiscal quantities
                # 1. total spin
                total_spin=sampling_magnetization(spin,N)
                mag=mag+total_spin/sweepnum
                mag2=mag2+(total_spin*total_spin)/sweepnum
                # 2. total energy
                total_energy=sampling_energy(spin,N,j0,ext)
                eng=eng+total_energy/sweepnum
                eng2=eng2+(total_energy*total_energy)/sweepnum
            else
            end
        end
        # cv, chi
        cv=(eng2-eng*eng)/(temp*temp)
        chi=(mag2-mag*mag)/(temp)
        return spin,mag,eng,cv,chi
    end

    # sampling energy 
    function sampling_energy(spin,N,j0,ext)
        total_energy=0
        for i in 1:N
            for j in 1:N
                if(i==N)
                    right=spin[1,j]
                else
                    right=spin[i+1,j]
                end
                if(j==N)
                    down=spin[i,1]
                else
                    down=spin[i,j+1]
                end
                total_energy=total_energy+(-j0)*spin[i,j]*(right+down)-ext*spin[i,j]
            end
        end
        return total_energy/(N*N)
    end

   # sampling magnetization
    function sampling_magnetization(spin,N)
        total_spin=0
        for i in 1:N
            for j in 1:N
                total_spin=total_spin+spin[i,j]
            end
        end
        return total_spin/(N*N)
    end

  # get energy
    function get_energy(spin,i,j,j0,N,ext)
        if(i==1)
            left  = spin[N,j]
            right = spin[2,j]
        elseif(i==N)
            left  = spin[N-1,j]
            right = spin[1,j]
        else
            left  = spin[i-1,j]
            right = spin[i+1,j]
        end
        if(j==1)
            up    = spin[i,N]
            down  = spin[i,2]
        elseif(j==N)
            up    = spin[i,N-1]
            down  = spin[i,1]
        else
            up    = spin[i,j-1]
            down  = spin[i,j+1]
        end
        
        deltaE = (-1.0)*j0*(spin[i,j]*(left+right+up+down))-spin[i,j]*ext
        return deltaE
    end
end #module
