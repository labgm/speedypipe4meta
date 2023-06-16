/* Courtesy of the answer provided by WhozCraig at StackOverflow
   https://stackoverflow.com/questions/23896421/efficiently-waiting-for-all-tasks-in-a-threadpool-to-finish
*/


#include "ThreadPool.h"

ThreadPool::ThreadPool(unsigned int n)	: busy(), processed(), stop()
//ThreadPool::ThreadPool(unsigned int n)
{
	for (unsigned int i=0; i<n; ++i)
		workers.emplace_back(std::bind(&ThreadPool::thread_proc, this));
}

ThreadPool::~ThreadPool()
{
	// set stop-condition
	std::unique_lock<std::mutex> latch(queue_mutex);
	stop = true;
	cv_task.notify_all();
	latch.unlock();

	// all threads terminate, then we're done.
	for (auto& t : workers)
		t.join();
}

void ThreadPool::thread_proc()
{
	while (true)
	{
		std::unique_lock<std::mutex> latch(queue_mutex);
		cv_task.wait(latch, [this](){ return stop || !tasks.empty(); });
		if (!tasks.empty())
		{
			// got work. set busy.
			++busy;

			// pull from queue
			auto fn = tasks.front();
			tasks.pop_front();

			// release lock. run async
			latch.unlock();

			// run function outside context
			fn();
			++processed;

			latch.lock();
			--busy;
			cv_finished.notify_one();
		}
		else if (stop)
		{
			break;
		}
	}
}

// generic function push
//template<class F> void ThreadPool::enqueue(F&& f)
void ThreadPool::enqueue(std::function <void (void)> func)
{
	std::unique_lock<std::mutex> lock(queue_mutex);
	tasks.emplace_back(std::forward<std::function<void(void)>>(func));
	cv_task.notify_one();
}

// waits until the queue is empty.
void ThreadPool::waitFinished()
{
	std::unique_lock<std::mutex> lock(queue_mutex);
	cv_finished.wait(lock, [this](){ return tasks.empty() && (busy == 0); });
}

