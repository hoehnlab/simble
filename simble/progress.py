import sys
import time

# def progressbar(it, prefix="", size=60, out=sys.stdout): # Python3.6+
#     count = len(it)
#     start = time.time() # time estimate start
#     def show(j):
#         x = int(size*j/count)
#         # time estimate calculation and string
#         remaining = ((time.time() - start) / j) * (count - j)        
#         mins, sec = divmod(remaining, 60) # limited to minutes
#         time_str = f"{int(mins):02}:{sec:03.1f}"
#         print(f"{prefix}[{u'█'*x}{('.'*(size-x))}] {j}/{count} Est wait {time_str}", end='\r', file=out, flush=True)
#     show(0.1) # avoid div/0 
#     for i, item in enumerate(it):
#         yield item
#         show(i+1)
#     print("\n", flush=True, file=out)

class ProgressBar():
    def __init__(
            self,
            count,
            prefix="",
            size=60,
            out=sys.stdout,
            line_to_overwrite=None,
    ):
        self.count = count
        self.size = size
        self.prefix = prefix
        self.out = out
        self.start = time.time()
        self.current = 0
        self.line_to_overwrite = line_to_overwrite
    
    def show(self, j):
        x = int(self.size*j/self.count)
        if j < 10:
            time_str = "Calculating wait..."
        else:
            remaining = ((time.time() - self.start) / j) * (self.count - j)
            mins, sec = divmod(remaining, 60)
            time_str = f"Estimated wait: {int(mins):02}:{sec:03.1f}"
        if self.line_to_overwrite:
            # read in everything in the output TextIOWrapper right now
            # and overwrite the line with the progress bar
            # lines = self.out.buffer.readlines()
            # lines[self.line_to_overwrite] = f"{self.prefix}[{u'█'*x}{('.'*(self.size-x))}] {j}/{self.count}  {time_str}\n"
            # # write the lines back to the output TextIOWrapper
            # self.out.writelines(lines)
            # print(' ' * (self.size+30), end='\r')
            print(f"\033[{self.line_to_overwrite};0H\033[K{self.prefix}[{u'█'*x}{('.'*(self.size-x))}] {j}/{self.count}  {time_str}", end='\n', file=self.out, flush=True)


            # with open(self.out, 'r') as f:
            #     lines = f.readlines()
            # lines[self.line_to_overwrite] = f"{self.prefix}[{u'█'*x}{('.'*(self.size-x))}] {j}/{self.count}  {time_str}\n"
            # with open(self.out, 'w') as f:
            #     f.writelines(lines)
        else:
            print(' ' * (self.size+30), end='\r')
            print(f"{self.prefix}[{u'█'*x}{('.'*(self.size-x))}] {j}/{self.count}  {time_str}", end='\r', file=self.out, flush=True)

    
    def update(self):
        self.current +=1
        self.show(self.current)
        if self.current >= self.count:
            print("\n", flush=True, file=self.out)