import java.util.*;


public class arrayFF implements java.util.Collection, java.io.Serializable, Cloneable, Indexed
    {
    public Object[] objs;
    public int numObjs;
    
    public arrayFF() { numObjs = 0; objs = new Object[1]; }
    
    /** Creates a DoublearrayFF with a given initial capacity. */
    public arrayFF(int capacity) { numObjs = 0; objs = new Object[capacity]; }
        
    /** Adds the objects from the other arrayFF without copying them.  The size of the
        new arrayFF is the minimum necessary size to hold the objects. */
    public arrayFF(final arrayFF other)
        {
        if (other==null) { numObjs = 0; objs = new Object[1]; }
        numObjs = other.numObjs;
        objs = new Object[numObjs];
        System.arraycopy(other.objs,0,objs,0,numObjs);
        }
    
    public int size()
        {
        return numObjs;
        }
    
    public boolean isEmpty()
        {
        return (numObjs<= 0);
        }
    
    public boolean addAll(final Collection other) { return addAll(numObjs, other.toArray()); }

    public boolean addAll(final int index, final Collection other) { return addAll(index, other.toArray()); }

    public boolean addAll(final int index, final Object[] other)
        {
        if (index < 0 || index > numObjs || 
            other==null || other.length == 0) return false;
        if (numObjs+other.length > objs.length)
            resize(2*(numObjs+other.length));
        if (index != numObjs)   
            System.arraycopy(objs,index,objs,index+other.length,other.length);
        System.arraycopy(other,0,objs,index,other.length);
        numObjs += other.length;
        return true;
        }
    
    public boolean addAll(final arrayFF other) { return addAll(numObjs,other); }

    public boolean addAll(final int index, final arrayFF other)
        {
        if (index < 0 || index > numObjs || 
            other==null || other.numObjs <= 0) return false;
        if (numObjs+other.numObjs > objs.length)
            resize(2*(numObjs+other.numObjs));
        if (index != numObjs)    
            System.arraycopy(objs,index,objs,index+other.numObjs,other.numObjs);
        System.arraycopy(other.objs,0,objs,index,other.numObjs);
        numObjs += other.numObjs;
        return true;
        }

    public Object clone() throws CloneNotSupportedException
        {
        arrayFF b = (arrayFF)(super.clone());
        b.objs = (Object[]) objs.clone();
        return b;
        }
    
    /** Resizes the internal array to at least the requested size. */
    public void resize(int toAtLeast)
        {
        if (objs.length >= toAtLeast)  
            return;

        if (objs.length * 2 > toAtLeast)  
            toAtLeast = objs.length * 2;

        // now resize
        Object[] newobjs = new Object[toAtLeast];
        System.arraycopy(objs,0,newobjs,0,numObjs);
        objs=newobjs;
        }
    
    /** Returns null if the arrayFF is empty, else returns the topmost object. */
    public Object top()
        {
        if (numObjs<= 0) return null;
        else return objs[numObjs-1];
        }
    
    /** Returns null if the arrayFF is empty, else removes and returns the topmost object. */
    public Object pop()
        {
        int numObjs = this.numObjs;
        if (numObjs<= 0) return null;
        Object ret = objs[--numObjs];
        objs[numObjs] = null; // let GC
        this.numObjs = numObjs;
        return ret;
        }
    

    public boolean push(final Object obj)
        {
        int numObjs = this.numObjs;
        if (numObjs >= objs.length) doubleCapacityPlusOne();
        objs[numObjs] = obj;
        this.numObjs = numObjs+1;
        return true;
        }
        
    public boolean add(final Object obj)
        {
        int numObjs = this.numObjs;
        if (numObjs >= objs.length) doubleCapacityPlusOne();
        objs[numObjs] = obj;
        this.numObjs = numObjs+1;
        return true;
        }
        

    void doubleCapacityPlusOne()
        {
        Object[] newobjs = new Object[numObjs*2+1];
        System.arraycopy(objs,0,newobjs,0,numObjs);
        objs=newobjs;
        }

    public boolean contains(final Object o)
        {
        int numObjs = this.numObjs;
        Object[] objs = this.objs;
        for(int x=0;x<numObjs;x++)
            if (o==null ?  objs[x]==null :  o==objs[x] || o.equals(objs[x])) return true;
        return false;
        }
        
    public boolean containsAll(final Collection c)
        {
        Iterator iterator = c.iterator();
        while(iterator.hasNext())
            if (!contains(iterator.next())) return false;
        return true;
        }

    public Object get(final int index)
        {
        if (index>=numObjs || index < 0)
            throwIndexOutOfBoundsException(index);
        return objs[index];
        }

    /** identical to get(index) */
    public Object getValue(final int index)
        {
        if (index>=numObjs || index < 0)
            throwIndexOutOfBoundsException(index);
        return objs[index];
        }

    public Object set(final int index, final Object element)
        {
        if (index>=numObjs || index < 0)
            throwIndexOutOfBoundsException(index);
        Object returnval = objs[index];
        objs[index] = element;
        return returnval;
        }

    /** identical to set(index, element) */
    public Object setValue(final int index, final Object element)
        {
        if (index>=numObjs || index < 0)
            throwIndexOutOfBoundsException(index);
        Object returnval = objs[index];
        objs[index] = element;
        return returnval;
        }

    public boolean removeAll(final Collection c)
        {
        boolean flag = false;
        Iterator iterator = c.iterator();
        while(iterator.hasNext())
            if (remove(iterator.next())) flag = true;
        return flag;
        }

    public boolean retainAll(final Collection c)
        {
        boolean flag = false;
        for(int x=0;x<numObjs;x++)
            if (!c.contains(objs[x]))
                {
                flag = true;
                remove(x);
                x--; 
                }
        return flag;
        }

    /** Removes the object at the given index, shifting the other objects down. */
    public Object removeNondestructively(final int index)
        {
        if (index>=numObjs || index < 0)
            throwIndexOutOfBoundsException(index);
        Object ret = objs[index];
        if (index < numObjs - 1)  
            System.arraycopy(objs, index+1, objs, index, numObjs - index - 1);
        objs[numObjs-1] = null;  
        numObjs--;
        return ret;
        }
    
    public boolean remove(final Object o)
        {
        for(int x=0;x<numObjs;x++)
            if (o==null ?  objs[x]==null :  o==objs[x] || o.equals(objs[x])) 
                {
                remove(x);
                return true;
                }
        return false;
        }
        
    /** Removes multiple instantiations of an object */
    public boolean removeMultiply(final Object o)
        {
        boolean flag = false;
        for(int x=0;x<numObjs;x++)
            if (o==null ?  objs[x]==null :  o==objs[x] || o.equals(objs[x])) 
                {
                flag = true;
                remove(x);
                x--;  
                }
        return flag;
        }

    /** Removes the object at the given index, moving the topmost object into its position. */
    public Object remove(final int index)
        {
        if (index>=numObjs || index < 0)
            throwIndexOutOfBoundsException(index);
        Object ret = objs[index];
        objs[index] = objs[numObjs-1];
        objs[numObjs-1] = null;  // let GC
        numObjs--;
        return ret;
        }
        
    protected void throwIndexOutOfBoundsException(final int index)
        {
        throw new IndexOutOfBoundsException(""+index);
        }
        
    public void clear()
        {
        numObjs = 0;
        }
        
    public Object[] toArray()
        {
        Object[] o = new Object[numObjs];
        System.arraycopy(objs,0,o,0,numObjs);
        return o;
        }
        
    public Object[] toArray(Object[] o)
        {
        if (o.length < numObjs)
            o = (Object[]) java.lang.reflect.Array.newInstance(o.getClass().getComponentType(), numObjs);
        System.arraycopy(objs,0,o,0,numObjs);
        if (o.length > numObjs)
            o[numObjs] = null;
        return null;
        }


    public Iterator iterator()
        {
        return new arrayFFIterator(this);
        }
    
    public Class componentType()
        {
        return null;
        }

    /** Sorts the arrayFF according to the provided comparator */
    public void sort(Comparator c) 
        {
        Arrays.sort(objs, 0, numObjs, c);
        }

    /** Replaces all elements in the arrayFF with the provided object. */
    public void fill(Object o)
        {
        // teeny bit faster
        Object[] objs = this.objs;
        int numObjs = this.numObjs;
        
        for(int x=0; x < numObjs; x++)
            objs[x] = o;
        }

    /** Shuffles (randomizes the order of) the arrayFF */
    public void shuffle(Random random)
        {
        // teeny bit faster
        Object[] objs = this.objs;
        int numObjs = this.numObjs;
        Object obj;
        int rand;
        
        for(int x=0; x < numObjs; x++)
            {
            rand = random.nextInt(numObjs);
            obj = objs[x];
            objs[x] = objs[rand];
            objs[rand] = obj;
            }
        }
    
    /** Shuffles (randomizes the order of) the arrayFF */
    public void shuffle(MersenneTwisterFast random)
        {
        // teeny bit faster
        Object[] objs = this.objs;
        int numObjs = this.numObjs;
        Object obj;
        int rand;
        
        for(int x=0; x < numObjs; x++)
            {
            rand = random.nextInt(numObjs);
            obj = objs[x];
            objs[x] = objs[rand];
            objs[rand] = obj;
            }
        }
    
    /** Reverses order of the elements in the arrayFF */
    public void reverse()
        {
        Object[] objs = this.objs;
        int numObjs = this.numObjs;
        int l = numObjs / 2;
        Object obj;
        for(int x=0; x < l; x++)
            {
            obj = objs[x];
            objs[x] = objs[numObjs - x - 1];
            objs[numObjs - x - 1] = obj;
            }
        }

    static class arrayFFIterator implements Iterator, java.io.Serializable
        {
        int obj = 0;
        arrayFF arrayFF;
        boolean canRemove = false;
        
        public arrayFFIterator(arrayFF arrayFF) { this.arrayFF = arrayFF; }
        
        public boolean hasNext()
            {
            return (obj < arrayFF.numObjs);
            }
        public Object next()
            {
            if (obj >= arrayFF.numObjs) throw new NoSuchElementException("No More Elements");
            canRemove = true;
            return arrayFF.objs[obj++];
            }
        public void remove()
            {
            if (!canRemove) throw new IllegalStateException("remove() before next(), or remove() called twice");
            if (obj - 1 >=  arrayFF.numObjs) throw new NoSuchElementException("No More Elements");
            arrayFF.removeNondestructively(obj-1);
            canRemove = false;
            }
        }
    }
