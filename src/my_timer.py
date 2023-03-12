import time

class My_Timer:
    def __init__(self):
        self.start()

    def start(self):
        self.start_time = time.perf_counter()
        self.last_time = self.start_time

    def print_now(self):
        print(f"{self.now():0.2f} seconds")

    def now(self):
        now_time = time.perf_counter()
        return now_time - self.start_time
    
    def point(self):
        """time from last call of point() function"""
        now_time = time.perf_counter()
        rez = now_time - self.last_time
        self.last_time = now_time
        return rez


def main():
    timer = My_Timer()

    rez = 0
    for i in range(3000):
        rez *= i ** i

    timer.print_now()


if __name__ == "__main__":
    main()