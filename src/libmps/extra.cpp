//  methods added to support the use of mpsolve by num programs

#include <mps/mps.h>
#include <Windows.h>

static bool bIsFatalExit;
static const char default_msg_folder[] = "\\temp";
static const char default_msg_path[] = "\\temp\\mpsolve_windows.txt";

#ifdef NUM_MESSAGING

#ifdef MPS_USE_PTHREADS
static mps_initialized_mutex(mpsolve_messaging_mutex_Static, PTHREAD_MUTEX_INITIALIZER);
#else
static std::mutex mpsolve_messaging_mutex_Static;
#endif
static bool mpsolve_messaging_folders_created;

//  make all parent folders for an output file

static void mps_make_parent_folders(const char* fullpath)
{
    int i;
    char ch;
    bool bHitFirst;
    char* mypath;

    mypath = (char*)mps_malloc(strlen(fullpath) + 1);
    strcpy(mypath, fullpath);

    bHitFirst = false;
    for (i = 0;;i++)
    {
        mypath[i] = toupper(mypath[i]);
        if (!mypath[i] || mypath[i] == '\\' || mypath[i] == '/')
        {
            if (bHitFirst)
            {
                ch = mypath[i];
                mypath[i] = 0;
                _mkdir(mypath);
                mypath[i] = ch;
            }
            if (!mypath[i]) break;
            bHitFirst = true;
        }
    }

    mps_free(mypath);
    mpsolve_messaging_folders_created = true;
}

// for debug purposes - output msg(s) to temp

int mps_msg_out(const char* msg, const char* msg_path)
{
    FILE* out;
    int i, tries;

    mps_mutex_lock(mpsolve_messaging_mutex_Static);

    if (!mpsolve_messaging_folders_created)
    {
        mps_make_parent_folders(default_msg_folder);
    }

    for (tries = 0;tries < 10;tries++)
    {
        out = fopen(msg_path, "a");
        if (out)
        {
            fseek(out, 0L, SEEK_END);
            for (i = 0;msg[i];i++);
            fwrite(msg, i, 1, out);
            fclose(out);
            mps_mutex_unlock(mpsolve_messaging_mutex_Static);
            return 0;
        }
        Sleep(10);
    }
    mps_mutex_unlock(mpsolve_messaging_mutex_Static);
    {
        fprintf(stderr, "Method Test_Issue_Msg_Global: failed to open msg file (%s)\n", strerror(errno));
        fprintf(stderr, "Method Test_Issue_Msg_Global: msg is \"%s\"\n", msg);
    }
    return(0);
}
#endif

//	   fatal error exit

void mps_fatal_exit(char* errmsg)
{
    char msg[MPS_ERRMSG_SIZE];

#ifndef NUM_USE_THREAD_SAFE
    while (bIsFatalExit)  // once a fatal exit occurs in mpsolve, pause all other mpsolve threads while this one exits
    {
        Sleep(1000);
    }
#endif

    bIsFatalExit = true;
    fprintf(stderr, "mpsolve_windows run-time error...\n");
    sprintf(msg, "%s\n", errmsg);
    fprintf(stderr, msg);
#ifdef NUM_MESSAGING
    mps_msg_out(msg, default_msg_path);
#endif

    exit(1);   // there are multiple threads, a 'throw' won't be be caught by any try/catch elsewhere
}

int mps_get_thread_id()
{
    return GetCurrentThreadId();
}

#ifndef MPS_USE_PTHREADS

static std::mutex guarded_list_mutex;

static struct MUTEX_OWNERSHIP
{
    struct MUTEX_OWNERSHIP* next_owner;
    mps_mutex_t* owned_mutex;
    int owner_thread_id;
} *first_owner;

//  lock but only if not already locked by the calling thread

bool mps_guarded_lock(mps_mutex_t& target_mutex)
{
    int sys_thread_id;
    MUTEX_OWNERSHIP* new_owner, * this_owner;
    bool bLockAcquired;

    sys_thread_id = mps_get_thread_id();
    bLockAcquired = false;

    while (1)
    {
        mps_mutex_lock(guarded_list_mutex);

        this_owner = first_owner;
        while (this_owner)
        {
            if (&target_mutex == this_owner->owned_mutex)
            {
                if (this_owner->owner_thread_id == sys_thread_id)
                {
                    mps_mutex_unlock(guarded_list_mutex);
                    if (bLockAcquired) return true;         // this method has just done the lock which was not a dup
                    return false;                           // the target mutex is already owned by this thread, a duplicate lock is not attempted
                }
                break;             // owned but by a different thread
            }
            this_owner = this_owner->next_owner;
        }
        if (!this_owner) // no one owns it, and I still own the list.  
        {
            // register myself as the owner so no one else can get this guarded lock
            // this will "lock" a guarded lock so long as user does not attempt a lock outside of the guard mechanism
            new_owner = new MUTEX_OWNERSHIP;
            new_owner->next_owner = first_owner;
            first_owner = new_owner;
            new_owner->owner_thread_id = sys_thread_id;
            new_owner->owned_mutex = &target_mutex;

            mps_mutex_unlock(guarded_list_mutex);

            // acquire the lock which probably won't wait since this thread is the registered owner
            mps_mutex_lock(target_mutex);
            bLockAcquired = true;
            continue;
        }

        // it is owned but by a different thread, ok to issue a lock which won't be a dup and will probably wait (after releasing the list)
        mps_mutex_unlock(guarded_list_mutex);
        mps_mutex_lock(target_mutex);   // use mutex logic to wait
        mps_mutex_unlock(target_mutex);
        continue;  // see if it is now available for locking after registering in the list
    }
    return false;  // won't get here
}

//  unlock but only if owned by the calling thread

bool mps_guarded_unlock(mps_mutex_t& target_mutex)
{
    int sys_thread_id;
    MUTEX_OWNERSHIP* last_owner, * this_owner;

    sys_thread_id = mps_get_thread_id();

    mps_mutex_lock(guarded_list_mutex);

    this_owner = first_owner;
    last_owner = NULL;
    while (this_owner)
    {
        if (&target_mutex == this_owner->owned_mutex)
        {
            if (this_owner->owner_thread_id == sys_thread_id)
            {
                // I own the guarded mutex, remove myself from the list
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
                mps_mutex_unlock(guarded_list_mutex);
                return true;         // this method has just done the release which was not a dup
            }
            break;             // owned but by a different thread
        }
        last_owner = this_owner;
        this_owner = this_owner->next_owner;
    }

    //  it was not owned by me, skip the release
    mps_mutex_unlock( guarded_list_mutex);
    return false;
}
#endif   // #ifndef MPS_USE_PTHREADS
