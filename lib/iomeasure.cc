#include "iomeasure.h"
#include <inttypes.h>
int IOMeasure::lru(const Array& a_array, const int a_cachesize)
{
    class carrier
    {
    public:
        const int   m_pageid;
        const int   m_accesstime;
    public:
        carrier(const int a_pid, const int a_time):
          m_pageid(a_pid), m_accesstime(a_time) {};
        ~carrier() {};
    };

    Hash cache;     // that keeps the access times of cached items
    Queue accessq;  // that order the accessed pages
    int pageread = 0;
    int timer = 1;
    for (int i=0; i<a_array.size(); i++)
    {
        int pid = (intptr_t)a_array.get(i);
        accessq.enqueue(new carrier(pid, timer));
        int latest = (intptr_t)cache.get(pid);
        if (latest == 0)
        {
            //-----------------------------------------------------------------
            // the page is admitted to the cache
            //-----------------------------------------------------------------
            cache.put(pid, (void*)timer);
            while (cache.size() > a_cachesize && !accessq.isEmpty())
            {
                //-------------------------------------------------------------
                // the cache is full, select victim to discard
                //-------------------------------------------------------------
                carrier* c = (carrier*)accessq.dequeue();
                latest = (intptr_t)cache.get(c->m_pageid);
                if (c->m_accesstime == latest)
                    cache.remove(c->m_pageid);
                delete c;
            }
            pageread++;
        }
        else
        {
            //-----------------------------------------------------------------
            // the page is already in the cache
            //-----------------------------------------------------------------
            cache.replace(pid, (void*)timer);
        }
        timer++;
    }

    //-------------------------------------------------------------------------
    // clean up
    //-------------------------------------------------------------------------
    while (!accessq.isEmpty())
        delete (carrier*)accessq.dequeue();
    cache.clean();

    return pageread;
}
