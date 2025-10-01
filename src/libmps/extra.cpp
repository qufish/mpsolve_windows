//  methods added to support the use of mpsolve by num programs

#include <mps/mps.h>
#include <Windows.h>

static bool bIsFatalExit;

//	   fatal error exit

void mps_fatal_exit(char* errmsg)
{
    char msg[ERRMSG_SIZE];

#ifndef NUM_USE_THREAD_SAFE
    while (bIsFatalExit)  // once a fatal exit occurs in mpsolve, pause all other mpsolve threads while this one exits
    {
        Sleep(1000);
    }
#endif

    bIsFatalExit = true;
    fprintf(stderr, "mpsolve_windows run-time error...\n");
    sprintf(msg, "%s\n", errmsg);
#ifdef NUM_MESSAGING
    issue_msg(msg);
#endif
    fprintf(stderr, msg);

    exit(1);   // there are multiple threads, a 'throw' won't be be caught by any try/catch elsewhere
}

static mps_mutex_t guard_list_mutex;
static bool bGuardListMutexInitialized;
static struct MUTEX_OWNERSHIP
{
    struct MUTEX_OWNERSHIP* next_owner;
    mps_mutex_t* owned_mutex;
    int owner_thread_id;
} *first_owner;

//  lock but only if not already locked by the calling thread

bool mps_tracked_lock(mps_mutex_t& target_mutex)
{
    int sys_thread_id;
    MUTEX_OWNERSHIP* new_owner, * this_owner;
    bool bLockAcquired;

    if (!bGuardListMutexInitialized)
    {
        mps_mutex_init(guard_list_mutex);
        bGuardListMutexInitialized = true;
    }

    sys_thread_id = GetCurrentThreadId();
    bLockAcquired = false;

    while (1)
    {
        mps_mutex_lock(guard_list_mutex);

        this_owner = first_owner;
        while (this_owner)
        {
            if (&target_mutex == this_owner->owned_mutex)
            {
                if (this_owner->owner_thread_id == sys_thread_id)
                {
                    mps_mutex_unlock(guard_list_mutex);
                    if (bLockAcquired) return true;         // this method has just done the lock which was not a dup
                    return false;                           // the target mutex is already owned by this thread, a duplicate lock is not attempted
                }
                break;             // owned but by a different thread
            }
            this_owner = this_owner->next_owner;
        }
        if (!this_owner) // no one owns it, and I still own the list.  
        {
            // register myself as the owner so no one else can get this tracked lock
            // this will "lock" a tracked lock so long as user does not attempt a lock outside of the guard mechanism
            new_owner = new MUTEX_OWNERSHIP;
            new_owner->next_owner = first_owner;
            first_owner = new_owner;
            new_owner->owner_thread_id = sys_thread_id;
            new_owner->owned_mutex = &target_mutex;

            mps_mutex_unlock(guard_list_mutex);

            // acquire the lock which probably won't wait since this thread is the registered owner
            mps_mutex_lock(target_mutex);
            bLockAcquired = true;
            continue;
        }

        // it is owned but by a different thread, ok to issue a lock which won't be a dup and will probably wait (after releasing the list)
        mps_mutex_unlock(guard_list_mutex);
        mps_mutex_lock(target_mutex);   // use mutex logic to wait
        mps_mutex_unlock(target_mutex);
        continue;  // see if it is now availble for locking after registering in the list
    }
    return false;  // won't get here
}

//  unlock but only if owned by the calling thread

bool mps_tracked_unlock(mps_mutex_t& target_mutex)
{
    int sys_thread_id;
    MUTEX_OWNERSHIP* last_owner, * this_owner;

    if (!bGuardListMutexInitialized)
    {
        mps_mutex_init(guard_list_mutex);
        bGuardListMutexInitialized = true;
    }

    sys_thread_id = GetCurrentThreadId();

    mps_mutex_lock(guard_list_mutex);

    this_owner = first_owner;
    last_owner = NULL;
    while (this_owner)
    {
        if (&target_mutex == this_owner->owned_mutex)
        {
            if (this_owner->owner_thread_id == sys_thread_id)
            {
                // I own the tracked mutex, remove myself from the list
                if (last_owner)
                {
                    last_owner->next_owner = this_owner->next_owner;
                }
                else
                {
                    first_owner = this_owner->next_owner;
                }
                delete this_owner;

                mps_mutex_unlock(target_mutex);
                mps_mutex_unlock(guard_list_mutex);
                return true;         // this method has just done the release which was not a dup
            }
            break;             // owned but by a different thread
        }
        last_owner = this_owner;
        this_owner = this_owner->next_owner;
    }

    //  it was not owned by me, skip the release
    mps_mutex_unlock( guard_list_mutex);
    return false;
}
