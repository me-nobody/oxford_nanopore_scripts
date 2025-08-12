import os
import socket
from mpi4py import MPI

user = os.environ['USER']
cpu = os.environ['BB_CPU']
world_comm = MPI.COMM_WORLD
world_size = world_comm.Get_size()
my_rank = world_comm.Get_rank()
node = socket.gethostname()

print(f"Hello {user}! I am process {my_rank} of {world_size} from {node} (node type: {cpu})")