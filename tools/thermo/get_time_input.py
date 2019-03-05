#Get Time Input Function

def get_time_input(time_start, time_final, time_bounce, N, N_bounce):
    import numpy as np
    time_bounce_floor = time_bounce - 0.5E-03;
    time_bounce_ceil  = time_bounce + 0.5E-03;

    N_prebounce  = 0.25 * N * ( time_bounce_floor - time_start ) / ( time_final - time_start ) 
    N_prebounce = int(N_prebounce)
    #N_prebounce = N_prebounce.astype(np.int32) doesn't work if not array
    N_postbounce = N - N_prebounce - N_bounce + 1

    time_input_bounce     = np.linspace(time_bounce_floor, time_bounce_ceil, N_bounce, endpoint=True)
    time_input_prebounce  = np.linspace(time_start, time_bounce_floor, N_prebounce, endpoint=True)
    time_input_postbounce = np.linspace(time_bounce_ceil, time_final, N_postbounce, endpoint=True)

    time_input = np.append(time_input_prebounce,time_bounce)
    time_input = np. concatenate((time_input, time_input_bounce, time_input_postbounce), axis=0)
    time_input = np.unique(time_input)
    
    return time_input