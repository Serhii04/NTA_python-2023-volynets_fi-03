import nta_algorithms as nalg

def main():
    while True:
        n = input("Type your number: ")
        if n.isnumeric():
            n = int(n)
            factor_n = nalg.get_canon_number_composition(n)

            if len(factor_n):
                print(f"n = {factor_n[0]}", end="")
                
                for d in range(1, len(factor_n)):
                    print(f" * {factor_n[d]}", end="")
                print("")
            
            return
        else:
            print(f"{n} isn't a number, please try again.")
            


if __name__ == "__main__":
    main()