from heapq import *
import random

# Flags to selectively enable parts of code
graphics = 0;
random_processing = 0;

if graphics:
  import numpy as np
  import matplotlib.pyplot as plt
  import matplotlib.colors as colors

# Constant definitions
EMPTY = 0;
FULL = 1;

FALSE = 0;
TRUE = 1;

# State definitions
IDLE              = 0;
INIT_SEND_LEFT    = 1;
INIT_SEND_RIGHT   = 2;
INIT_SEND_UP      = 3;
INIT_SEND_DOWN    = 4;
SEND_LEFT         = 5;
SEND_RIGHT        = 6;
SEND_UP           = 7;
SEND_DOWN         = 8;
RECV_LEFT         = 9;
RECV_RIGHT        = 10;
RECV_UP           = 11;
RECV_DOWN         = 12;

# Event definitions
START_LEFT_RIGHT  = 1;
START_UP_DOWN     = 2;
SEND_ACK_ARRIVE   = 3;
SEND_NACK_ARRIVE  = 4;
SEND_COMPLETE     = 5;
PKT_HEAD_ARRIVE   = 6;
PKT_TAIL_ARRIVE   = 7;

# Event source definitions
LEFT              = 0;
SRC_LEFT          = 0;
RIGHT             = 1;
SRC_RIGHT         = 1;
UP                = 2;
SRC_UP            = 2;
DOWN              = 3;
SRC_DOWN          = 3;
SRC_SELF          = 4;

# Network parameters
Nx = 64;
Ny = 64;
Nf = 256;
pkt_size = Nf * 2 * 4;  # in bytes
link_speed = 1e9;       # bits per second
latency = 4e-9;
bw_delay = (pkt_size * 8) / link_speed;
bw_delay_1col = (pkt_size * 128 * 8) / link_speed;
resend_delay = bw_delay / 2;
NUM_CYCLES = 3;

# Random processing parameters
mean_time = 2e-6;
std_dev = 0.5e-6;
min_time = 1e-6;
max_time = 4e-6;

# Define the event class
# NOTE: the source is defined as the ORIGINAL source of the original message.
class Event:

  def __init__(self, event_type, event_time, source):
    self.event_time = event_time;
    self.event_type = event_type;
    self.source = source;

  def __cmp__(self, other):
    if self.event_time < other.event_time:
      return -1;
    elif self.event_time == other.event_time:
      return 0;
    else:
      return 1;


# Function to print the given state
def print_state(state):
  if state == IDLE: return 'IDLE';
  elif state == INIT_SEND_LEFT: return 'INIT_SEND_LEFT';
  elif state == INIT_SEND_RIGHT: return 'INIT_SEND_RIGHT';
  elif state == INIT_SEND_UP: return 'INIT_SEND_UP';
  elif state == INIT_SEND_DOWN: return 'INIT_SEND_DOWN';
  elif state == SEND_LEFT: return 'SEND_LEFT';
  elif state == SEND_RIGHT: return 'SEND_RIGHT';
  elif state == SEND_UP: return 'SEND_UP';
  elif state == SEND_DOWN: return 'SEND_DOWN';
  elif state == RECV_LEFT: return 'RECV_LEFT';
  elif state == RECV_RIGHT: return 'RECV_RIGHT';
  elif state == RECV_UP: return 'RECV_UP';
  elif state == RECV_DOWN: return 'RECV_DOWN';

# Function to print the type and source of the given event
def print_event(e):
  ret_str = '';
  if e.event_type == START_LEFT_RIGHT: ret_str = 'START_LEFT_RIGHT';
  elif e.event_type == START_UP_DOWN: ret_str = 'START_UP_DOWN';
  elif e.event_type == SEND_ACK_ARRIVE: ret_str = 'SEND_ACK_ARRIVE';
  elif e.event_type == SEND_NACK_ARRIVE: ret_str = 'SEND_NACK_ARRIVE';
  elif e.event_type == SEND_COMPLETE: ret_str = 'SEND_COMPLETE';
  elif e.event_type == PKT_HEAD_ARRIVE: ret_str = 'PKT_HEAD_ARRIVE';
  elif e.event_type == PKT_TAIL_ARRIVE: ret_str = 'PKT_TAIL_ARRIVE';

  if e.source == LEFT: ret_str += ' from LEFT';
  elif e.source == RIGHT: ret_str += ' from RIGHT';
  elif e.source == UP: ret_str += ' from UP';
  elif e.source == DOWN: ret_str += ' from DOWN';
  elif e.source == SRC_SELF: ret_str += ' from SELF';

  return ret_str;

# Function to generate random time
def rand_time():
  t = random.gauss(mean_time, std_dev);
  if t < min_time: t = min_time;
  if t > max_time: t = max_time;
  return t;


# Define each node in the network
class Node:

  # Constructor
  def __init__(self, startup_delay, ant_x, ant_y):
    self.state = IDLE;
    self.ant_x = ant_x;
    self.ant_y = ant_y;
    self.event_queue = [];

    self.send_count = 0;
    self.send_complete_count = 0;
    self.dir_msg_count = [0, 0, 0, 0];
    self.cycle = 0;
    self.finished = FALSE;
    self.send_buf = [EMPTY, EMPTY, EMPTY, EMPTY];
    self.recv_buf = [EMPTY, EMPTY, EMPTY, EMPTY];
    self.recv_pending = [FALSE, FALSE, FALSE, FALSE];
    self.recv_complete = [FALSE, FALSE, FALSE, FALSE];

    self.add_event( Event(START_LEFT_RIGHT, startup_delay, SRC_SELF) );

  # Add an event to the node
  def add_event(self, e):
    heappush(self.event_queue, e);

  # Get the time of the next event
  def next_event_time(self):
    if len(self.event_queue) > 0:
      return self.event_queue[0].event_time;
    else:
      return -1;

  # Determine the opposite direction to the provided direction
  def opp(self, direction):
    if   direction == LEFT:   return RIGHT;
    elif direction == RIGHT:  return LEFT;
    elif direction == UP:     return DOWN;
    elif direction == DOWN:   return UP;
    else:                   
      print 'ERROR: no opposite direction.';
      return -1;

  # Create events and update buffers for sending a message
  def send_msg(self, direction, neighbors):
    neighbors[direction].add_event( Event(PKT_HEAD_ARRIVE, time + latency, self.opp(direction)) );
    self.send_count += 1;
    self.send_buf[direction] = FULL;

  # Determines the bandwidth delay, based on the current cycle
  def get_bw_delay(self):
    return bw_delay if self.cycle == 0 else bw_delay_1col;

  # Fire the next event
  def fire(self, time, neighbors):

    next_event = heappop(self.event_queue);

    if self.ant_x == 0 and self.ant_y == 2:
      print 'Time:', time, ' Node [0][2]: state =', print_state(self.state), ' event =', print_event(next_event), \
            ' msg_count =', sum(self.dir_msg_count), ' cycle=', self.cycle, ' lcount=', self.dir_msg_count[LEFT], ' rcount=', self.dir_msg_count[RIGHT];

    # Start event for L-R communication
    if next_event.event_type == START_LEFT_RIGHT:
      if neighbors[LEFT] != -1:
        self.send_msg(LEFT, neighbors);
        self.state = INIT_SEND_LEFT;
      else:
        self.state = RECV_RIGHT;

    # Start event for U-D communication
    elif next_event.event_type == START_UP_DOWN:
      if neighbors[UP] != -1:
        self.send_msg(UP, neighbors);
        self.state = INIT_SEND_UP;
      else:
        self.send_msg(DOWN, neighbors);
        self.state = INIT_SEND_DOWN;

    # Event: acknowledgement of send is received
    elif next_event.event_type == SEND_ACK_ARRIVE:

      if self.state == INIT_SEND_LEFT:
          self.state = RECV_RIGHT;

      if self.state == INIT_SEND_UP:
        if neighbors[DOWN] != -1:
          self.send_msg(DOWN, neighbors);
          self.state = INIT_SEND_DOWN;    
        else:
          self.state = RECV_UP;

      elif self.state == INIT_SEND_DOWN:
        self.state = RECV_DOWN;

      elif self.state == SEND_LEFT:
          self.state = RECV_RIGHT;

      elif self.state == SEND_RIGHT:
        if self.dir_msg_count[RIGHT] < (Nx-1 - self.ant_x):
          self.state = RECV_RIGHT;
        else:
          self.state = RECV_LEFT;

      elif self.state == SEND_UP:
        if self.dir_msg_count[UP] < self.ant_y:
          self.state = RECV_UP;
        else:
          self.state = RECV_DOWN;

      elif self.state == SEND_DOWN:
        if self.dir_msg_count[DOWN] < (Ny-1 - self.ant_y):
          self.state = RECV_DOWN;
        else:
          self.state = RECV_UP;

    # Event: SEND finishes, and so send buffer is now available
    elif next_event.event_type == SEND_COMPLETE:
      self.send_buf[self.opp(next_event.source)] = EMPTY;
      self.send_complete_count += 1;

    # Event: head of a data packet is received
    elif next_event.event_type == PKT_HEAD_ARRIVE:

      # If space in buffer
      if self.recv_buf[next_event.source] == EMPTY:
        self.add_event( Event(PKT_TAIL_ARRIVE, time + self.get_bw_delay(), next_event.source) );
        neighbors[next_event.source].add_event( Event(SEND_ACK_ARRIVE, time + latency, next_event.source) );
        neighbors[next_event.source].add_event( Event(SEND_COMPLETE, time - latency + self.get_bw_delay(), next_event.source) );
        self.recv_buf[next_event.source] = FULL;

      # Else no space in buffer
      else:
        self.recv_pending[next_event.source] = TRUE;

    # Event: tail of a data packet is received
    elif next_event.event_type == PKT_TAIL_ARRIVE:
      self.recv_complete[next_event.source] = TRUE;

    # If we are blocked on a receive, check if message is fully received.
    for direction in range(4):

      if self.state == (RECV_LEFT + direction) and self.recv_complete[direction]:

        # If we can receive the message
        # That is, If no neighbor on the opposite 'direction' side or
        # If space in the opposite send buffer
        if neighbors[self.opp(direction)] == -1 or self.send_buf[self.opp(direction)] == EMPTY:
          self.recv_buf[direction] = EMPTY;
          self.recv_complete[direction] = FALSE;
          self.dir_msg_count[direction] += 1;

          # If we have a neighbor in the opposite direction, forward message on
          if neighbors[self.opp(direction)] != -1:
            self.send_msg(self.opp(direction), neighbors);
            self.state = (SEND_LEFT + self.opp(direction));

          # If a receive was pending (because receive buffer was full), send an ACK and start receiving next msg
          if self.recv_pending[direction]: 
            self.add_event( Event(PKT_TAIL_ARRIVE, time + 2*latency + self.get_bw_delay(), direction) );
            neighbors[direction].add_event( Event(SEND_ACK_ARRIVE, time + latency, direction) );
            neighbors[direction].add_event( Event(SEND_COMPLETE, time + latency + self.get_bw_delay(), direction) );
            self.recv_pending[direction] = FALSE;
            self.recv_buf[direction] = FULL;

        # Else no space in the 'opposite' buffer
        else:
          pass;
          #print 'ERROR: cannot forward message on due to insufficient space in own send buffer', self.ant_x, self.ant_y, self.cycle;

    # If we have received all the messages, finish    
    if ( (self.cycle == 0 and self.dir_msg_count[RIGHT] == (Nx-1 - self.ant_x)) or       \
         (self.cycle >= 1 and self.dir_msg_count[UP] + self.dir_msg_count[DOWN] == (Ny - 1)) ) and \
       self.send_count == self.send_complete_count:

      self.cycle += 1;

      # If we have completed all 3 cycles of communication
      if self.cycle >= NUM_CYCLES or self.ant_x > 0:
        self.finished = TRUE;

      # Else do some computation, and then move to next phase
      else:
        # Reset counters
        self.send_count = 0;
        self.send_complete_count = 0;
        if self.cycle == 1:
          self.dir_msg_count[LEFT] = 0;
          self.dir_msg_count[RIGHT] = 0; 
        else:
          self.dir_msg_count[UP] = 0;
          self.dir_msg_count[DOWN] = 0; 

        # Start next cycle after computation is complete
        if random_processing:
          computation_delay = rand_time();
        else:
          computation_delay = mean_time;

        self.add_event( Event(START_UP_DOWN, time + computation_delay, SRC_SELF) );

# Create the network
network = []
for i in range(Nx):
  net = [];

  for j in range(Ny):
    net.append( Node(startup_delay = rand_time() if random_processing else mean_time, \
                    ant_x = i, ant_y = j) );

  network.append(net);

time = 0;
last_plot_time = 0;
global_next_event = 0;
done = FALSE;

# Setup the graphics simulation
if graphics:
  plt.ion()
  normalizer = colors.Normalize(vmin=0, vmax=NUM_CYCLES, clip=False);
  cycle_data = [[node.cycle for node in row ] for row in network]
  mat_data = np.array(cycle_data)
  plt.pcolor(mat_data, norm=normalizer)
  plt.draw()

# Step through time until simulation finished
while(not done):
  
  done = TRUE;

  if global_next_event != 9999:
    time = global_next_event;

  global_next_event = 9999;

  # Iterate through all nodes
  for x in range(Nx):
    for y in range(Ny):

      node = network[x][y];

      # If next_event_time == time, fire the node
      if node.next_event_time() == time:
        left_n  = network[x-1][y] if x > 0      else -1;
        right_n = network[x+1][y] if x < Nx - 1 else -1;
        up_n    = network[x][y-1] if y > 0      else -1;
        down_n  = network[x][y+1] if y < Ny - 1 else -1;

        neighbors = [left_n, right_n, up_n, down_n];
        node.fire(time, neighbors);

      # If node not finished, make done FALSE
      if not node.finished:
        done = FALSE;

      # If the next_event_time is less than current global_next_event, update global_next_event
      if node.next_event_time() < global_next_event and node.next_event_time() != -1:
        global_next_event = node.next_event_time();

  # Update the plot of the array
  if graphics and (time - last_plot_time) > 10e-3:
    cycle_data = [[node.cycle for node in row ] for row in network]
    mat_data = np.array(cycle_data)
    plt.pcolor(mat_data.T, norm=normalizer)
    plt.draw()
    last_plot_time = time;

theoretical_time = NUM_CYCLES * mean_time + (bw_delay + latency) * (Nx-1) + \
                    (NUM_CYCLES-1) * (bw_delay_1col + latency) * (Ny-1);
print 'Simulation finished at time:     ', time, 'seconds.';
print 'Theoretical completion time is:  ', theoretical_time, 'seconds';

if graphics:
  plt.ioff()
  plt.show()


