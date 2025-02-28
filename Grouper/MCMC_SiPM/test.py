import subprocess
from multiprocessing import shared_memory
import struct

# Run the C++ program

#hide cout of the program
subprocess.run("./MCMC_SiPM 50 30 0.68 4.0 1e-5 1", shell=True, stdout=subprocess.PIPE)

# Attach to the shared memory
shm_name = "1"
shm = shared_memory.SharedMemory(name="1")

# Read the double value
data = shm.buf[:8]
value = struct.unpack("d", data)[0]

print(f"Value read from shared memory: {value}")


