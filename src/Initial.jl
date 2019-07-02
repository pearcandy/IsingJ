module Initial

    export random_start,cold_start

    #random start : initialization of spin
    function random_start(N)
        spin=ones(N,N)
        for i in 1:N
            for j in 1:N
                if(rand()<=0.5)
                    spin[i,j]=-1
                else
                    spin[i,j]=1
                end
            end
        end        
        return spin
    end
 
    #cold start : initialization of spin
    function cold_start(N)
        spin=ones(N,N)
        return spin
    end
end # module
