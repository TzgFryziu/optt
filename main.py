from opt_alg import *
from user_fun import *
import numpy as np
import csv

# def lab1():
#
#     result = {"exp":[],"fib":[],"lag":[]}
#
#     for i in range(100):
#         a,b,x0,fcalls = expansion(ff1T, np.random.uniform(-100,100), 0.1, 3, 100)
#         result["exp"].append([a,b,x0,fcalls])
#         result["fib"].append(fib(ff1T, -100,100, 0.001))
#         result["lag"].append(lag(ff1T, -100, 100, np.random.uniform(a,b) , 0.001, 0.001, 100))
#
#     # Zapis do pliku CSV
#     with open("wynik_bez_zawezenia.csv", "w", newline='') as csv_file:
#         writer = csv.writer(csv_file)
#         # Zapisujemy nagłówki
#         writer.writerow(["Index","Exp x0","Exp a","Exp b","Exp fcalls", "Fib x","Fib y","Fib fcalls", "Lag x", "Lag y","Lag fcalls"])
#
#         # Zapisujemy dane
#         for i in range(100):
#             writer.writerow([i,result["exp"][i][2],result["exp"][i][0],result["exp"][i][1],result["exp"][i][3], result["fib"][i][0],result["fib"][i][1],result["fib"][i][2], result["lag"][i][0],result["lag"][i][1],result["lag"][i][2]])
#
def lab1():
    # Define parameters for the tank system
    P_A = 0.5
    V_0A = 5
    T_0A = 90
    P_B = 1
    V_0B = 1
    T_0B = 20
    F_in_B = 0.01
    T_in_B = 20
    D_B = 0.00365665
    a = 0.98
    b = 0.63
    g = 9.81

    # Parameter list for the differential equation function
    params = [P_A,V_0A,T_0A,P_B,V_0B,T_0B,F_in_B,T_in_B,D_B,a,b,g]


    D_A_interval = [0.0001,0.01]

    print("Interval for D_A after expansion search:", D_A_interval)

    # Step 2: Apply Fibonacci method within the interval
    epsilon = 0.00001
    gamma = 0.00001
    # optimal_D_A = fib_r(max_temp_in_tank, D_A_interval[0], D_A_interval[1], epsilon,params)

    # print("Optimal D_A:", optimal_D_A)
    # # optimal_D_A = 0.001166
    # max_temp = max_temp_in_tank(optimal_D_A, params)
    # print(f"Maximum Temperature in Tank B at optimal D_A: {max_temp+50} °C")
    #
    optimal_D_A = lag_r(max_temp_in_tank,D_A_interval[0],D_A_interval[1],epsilon,gamma,1000,params)
    # print("Optimal D_A:", optimal_D_A)
    # # optimal_D_A = 0.001166
    # max_temp = max_temp_in_tank(optimal_D_A, params)
    # print(f"Maximum Temperature in Tank B at optimal D_A: {max_temp+50} °C")
    #



if __name__ == '__main__':
    
    lab1()
    

