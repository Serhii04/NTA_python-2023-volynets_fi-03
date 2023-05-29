import nta_algorithms_lab_3 as nalg


def main():

    print("To solve a^x = b (mod p), type next coeficients:")


    alpha = int(input("Type a value: a = "))
    beta = int(input("Type b value: b = "))
    p = int(input("Type p value: p = "))

    nalg.index_calculus_timed(alpha=alpha, beta=beta, n=p - 1)
    

            


if __name__ == "__main__":
    main()