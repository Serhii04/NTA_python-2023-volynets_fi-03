import nta_algorithms_lab_2 as nalg
from nta_algorithms_lab_2 import time_limit
import my_timer
import signal

def main():

    print("To solve a^x = b (mod p), type next coeficients:")

    timeout = 60*5

    alpha = int(input("Type a value: a = "))
    beta = int(input("Type b value: b = "))
    p = int(input("Type p value: p = "))

    rez = None
    timer = my_timer.My_Timer()
    try:
        with time_limit(timeout):
            rez = nalg.SPG(alpha=alpha, beta=beta, n=p - 1)
    except UserWarning as e:
        print("Timed out!")
    
    if rez:
        # print(f"{alpha}^x = {beta} (mod {n+1})")
        print(f"x = {rez}, spend: {timer.now():0.4f}s")
    else:
        print(f"No result, spend: {timer.now():0.4f}s")

    return rez
            


if __name__ == "__main__":
    main()