# University Project - Parallel Merge Sort with OpenMP

## Overview

This project implements a **parallel Merge Sort algorithm** using **OpenMP**. The main objectives are to:
- **Follow the code hints:** Utilize provided hints to guide the implementation.
- **Use Tasks with Synchronization:** Apply OpenMP tasks when needed and add appropriate synchronization.
- **Create a Parallel Region:** Establish a parallel region with a single creator and multiple executors.
- **Implement a Cut-off Mechanism:** Integrate a cut-off threshold to switch from parallel to sequential execution when advantageous.
- **Test Different Configurations:** Experiment with various configurations to assess performance improvements.

## Project Requirements

- **Parallel Merge Sort Implementation:** Develop the merge sort algorithm leveraging OpenMP for concurrent execution.
- **Task-Based Parallelism:** Use OpenMP tasks for recursive sorting and ensure proper synchronization among tasks.
- **Single Creator with Multiple Executors:** Design the parallel region such that a single thread creates tasks while multiple threads execute them.
- **Cut-off Mechanism:** Add a mechanism to determine when to switch to a sequential sort to minimize overhead.
- **Configuration Testing:** Provide methods to test various configurations (e.g., different cut-off thresholds and thread counts) to optimize performance.


