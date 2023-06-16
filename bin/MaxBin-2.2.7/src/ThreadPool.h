#ifndef __THREADPOOL_H__
#define __THREADPOOL_H__

/* Courtesy of the answer provided by WhozCraig at StackOverflow
   https://stackoverflow.com/questions/23896421/efficiently-waiting-for-all-tasks-in-a-threadpool-to-finish
*/

#include "global_inc.h"
#include <atomic>
#include <vector>
#include <deque>
#include <functional>
#include <thread>
#include <condition_variable>
#include <mutex>

//thread pool
class ThreadPool
{
	public:
		ThreadPool(unsigned int n = std::thread::hardware_concurrency());

		//template<class F> void enqueue(F&& f);
		void enqueue(std::function <void (void)> func);
		void waitFinished();
		~ThreadPool();

		unsigned int getProcessed() const { return processed; }

	private:
		void thread_proc();

		std::vector< std::thread > workers;
		std::deque< std::function<void()> > tasks;
		std::mutex queue_mutex;
		std::condition_variable cv_task;
		std::condition_variable cv_finished;
		std::atomic_uint processed;
		unsigned int busy;
		bool stop;

};

#endif

