import java.io.*;
import java.util.*;

/** 
* MersenneTwisterFast is a drop-in subclass replacement for java.util.Random. It is properly synchronized and can be used in a multithreaded environment.
* MersenneTwisterFast is not a subclass of java.util.Random. It has the same public methods as Random does, however, and it is algorithmically identical 
* to MersenneTwister. MersenneTwisterFast has hard-code inlined all of its methods directly, and made all of them final (well, the ones of consequence anyway). 
* Further, these methods are not synchronized, so the same MersenneTwisterFast instance cannot be shared by multiple threads. But all this helps 
* MersenneTwisterFast achieve over twice the speed of MersenneTwister.
* About the Mersenne Twister. This is a Java version of the C-program for MT19937: Integer version. next(32) generates one pseudorandom unsigned integer 
* (32bit) which is uniformly distributed among 0 to 2^32-1 for each call. next(int bits) >>>'s by (32-bits) to get a value ranging between 0 and 2^bits-1 
* long inclusive; hope that's correct. setSeed(seed) set initial values to the working area of 624 words. For setSeed(seed), seed is any 32-bit integer except for 0.
* Orignally Coded by Takuji Nishimura, considering the suggestions by Topher Cooper and Marc Rieffel in July-Aug. 1997
* Translated to Java by Michael Lecuyer January 30, 1999 Copyright (C) 1999 Michael Lecuyer
* This library is free software; you can redistribute it and or modify it under the terms of the GNU Library General Public License as published by the Free 
* Software Foundation; either version 2 of the License, or (at your option) any later version. This library is distributed in the hope that it will be useful, 
* but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Library General Public License 
* for more details. You should have received a copy of the GNU Library General Public License along with this library; if not, write to the Free Foundation, Inc., 
* 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
* Makoto Matsumoto and Takuji Nishimura, the original authors ask "When you use this, send an email to: matumoto@math.keio.ac.jp with an appropriate reference to 
* your work" You might also point out this was a translation.
* Reference. M. Matsumoto and T. Nishimura, "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform Pseudo-Random Number Generator", ACM Transactions on 
* Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, pp 3--30.
* Bug Fixes. This implementation implements the bug fixes made in Java 1.2's version of Random, which means it can be used with earlier versions of Java. See the 
* JDK 1.2 java.util.Random documentation for further documentation on the random-number generation contracts made. Additionally, there's an undocumented bug in 
* the JDK java.util.Random.nextBytes() method, which this code fixes.
* Important Note. Just like java.util.Random, this generator accepts a long seed but doesn't use all of it. java.util.Random uses 48 bits. The Mersenne Twister 
* instead uses 32 bits (int size). So it's best if your seed does not exceed the int range.
**/

public class MersenneTwisterFast implements Serializable
    {
    // Period parameters
    private static final int N = 624;
    private static final int M = 397;
    private static final int MATRIX_A = 0x9908b0df;   
    private static final int UPPER_MASK = 0x80000000; 
    private static final int LOWER_MASK = 0x7fffffff; 


    // Tempering parameters
    private static final int TEMPERING_MASK_B = 0x9d2c5680;
    private static final int TEMPERING_MASK_C = 0xefc60000;
    
    private int mt[]; 
    private int mti; 
    private int mag01[];
    
    //private static final long GOOD_SEED = 4357;

    private double __nextNextGaussian;
    private boolean __haveNextNextGaussian;


    /**
     * Constructor using the default seed.
     */
    public MersenneTwisterFast()
        {
	this(System.currentTimeMillis());
        }
    
    /**
     * Constructor using a given seed.  
     *
     */
    public MersenneTwisterFast(final long seed)
        {
	setSeed(seed);
        }
    

    /**
     * Constructor using an array.
     */
    public MersenneTwisterFast(final int[] array)
        {
	setSeed(array);
        }


    /**
     * Initialize the pseudo random number generator.  
     */

    synchronized public void setSeed(final long seed)
        {

	__haveNextNextGaussian = false;

	mt = new int[N];
	
	mag01 = new int[2];
	mag01[0] = 0x0;
	mag01[1] = MATRIX_A;

        mt[0]= (int)(seed & 0xfffffff);
        for (mti=1; mti<N; mti++) 
            {
            mt[mti] = 
		(1812433253 * (mt[mti-1] ^ (mt[mti-1] >>> 30)) + mti); 

            mt[mti] &= 0xffffffff;
            }
        }


    /**
     * An alternative, more complete, method of seeding the
     * pseudo random number generator.  
     */

    synchronized public void setSeed(final int[] array)
	{
        int i, j, k;
        setSeed(19650218);
        i=1; j=0;
        k = (N>array.length ? N : array.length);
        for (; k!=0; k--) 
            {
            mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >>> 30)) * 1664525)) + array[j] + j; /* non linear */
            mt[i] &= 0xffffffff; 
            i++;
            j++;
            if (i>=N) { mt[0] = mt[N-1]; i=1; }
            if (j>=array.length) j=0;
            }
        for (k=N-1; k!=0; k--) 
            {
            mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >>> 30)) * 1566083941)) - i; /* non linear */
            mt[i] &= 0xffffffff; 
            i++;
            if (i>=N) 
                {
                mt[0] = mt[N-1]; i=1; 
                }
            }
        mt[0] = 0x80000000; 
        }


    public final int nextInt()
	{
	int y;
	
	if (mti >= N)   // generate N words at one time
	    {
	    int kk;
            final int[] mt = this.mt; // locals are slightly faster 
            final int[] mag01 = this.mag01; // locals are slightly faster 
	    
	    for (kk = 0; kk < N - M; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+M] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    for (; kk < N-1; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+(M-N)] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
	    mt[N-1] = mt[M-1] ^ (y >>> 1) ^ mag01[y & 0x1];

	    mti = 0;
	    }
  
	y = mt[mti++];
	y ^= y >>> 11;                          // TEMPERING_SHIFT_U(y)
	y ^= (y << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(y)
	y ^= (y << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(y)
	y ^= (y >>> 18);                        // TEMPERING_SHIFT_L(y)

	return y;
	}



    public final short nextShort()
	{
	int y;
	
	if (mti >= N)   // generate N words at one time
	    {
	    int kk;
            final int[] mt = this.mt; // locals are slightly faster 
            final int[] mag01 = this.mag01; // locals are slightly faster 
	    
	    for (kk = 0; kk < N - M; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+M] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    for (; kk < N-1; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+(M-N)] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
	    mt[N-1] = mt[M-1] ^ (y >>> 1) ^ mag01[y & 0x1];

	    mti = 0;
	    }
  
	y = mt[mti++];
	y ^= y >>> 11;                          // TEMPERING_SHIFT_U(y)
	y ^= (y << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(y)
	y ^= (y << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(y)
	y ^= (y >>> 18);                        // TEMPERING_SHIFT_L(y)

	return (short)(y >>> 16);
	}



    public final char nextChar()
	{
	int y;
	
	if (mti >= N)   // generate N words at one time
	    {
	    int kk;
            final int[] mt = this.mt; // locals are slightly faster 
            final int[] mag01 = this.mag01; // locals are slightly faster 
	    
	    for (kk = 0; kk < N - M; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+M] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    for (; kk < N-1; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+(M-N)] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
	    mt[N-1] = mt[M-1] ^ (y >>> 1) ^ mag01[y & 0x1];

	    mti = 0;
	    }
  
	y = mt[mti++];
	y ^= y >>> 11;                          // TEMPERING_SHIFT_U(y)
	y ^= (y << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(y)
	y ^= (y << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(y)
	y ^= (y >>> 18);                        // TEMPERING_SHIFT_L(y)

	return (char)(y >>> 16);
	}


    public final boolean nextBoolean()
	{
	int y;
	
	if (mti >= N)   // generate N words at one time
	    {
	    int kk;
            final int[] mt = this.mt; // locals are slightly faster 
            final int[] mag01 = this.mag01; // locals are slightly faster 
	    
	    for (kk = 0; kk < N - M; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+M] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    for (; kk < N-1; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+(M-N)] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
	    mt[N-1] = mt[M-1] ^ (y >>> 1) ^ mag01[y & 0x1];

	    mti = 0;
	    }
  
	y = mt[mti++];
	y ^= y >>> 11;                          // TEMPERING_SHIFT_U(y)
	y ^= (y << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(y)
	y ^= (y << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(y)
	y ^= (y >>> 18);                        // TEMPERING_SHIFT_L(y)

	return (boolean)((y >>> 31) != 0);
	}



    public final boolean nextBoolean(final float probability)
	{
	int y;
	
	if (probability < 0.0f || probability > 1.0f)
	    throw new IllegalArgumentException ("probability must be between 0.0 and 1.0 inclusive.");
	if (probability==0.0f) return false;		
	else if (probability==1.0f) return true;	
	if (mti >= N)   
	    {
	    int kk;
            final int[] mt = this.mt; 
            final int[] mag01 = this.mag01; 
	    
	    for (kk = 0; kk < N - M; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+M] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    for (; kk < N-1; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+(M-N)] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
	    mt[N-1] = mt[M-1] ^ (y >>> 1) ^ mag01[y & 0x1];

	    mti = 0;
	    }
  
	y = mt[mti++];
	y ^= y >>> 11;                          // TEMPERING_SHIFT_U(y)
	y ^= (y << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(y)
	y ^= (y << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(y)
	y ^= (y >>> 18);                        // TEMPERING_SHIFT_L(y)

	return (y >>> 8) / ((float)(1 << 24)) < probability;
	}



    public final boolean nextBoolean(final double probability)
	{
	int y;
	int z;

	if (probability < 0.0 || probability > 1.0)
	    throw new IllegalArgumentException ("probability must be between 0.0 and 1.0 inclusive.");
	if (probability==0.0) return false;		
	else if (probability==1.0) return true;	
	if (mti >= N)   
	    {
	    int kk;
            final int[] mt = this.mt; // locals are slightly faster 
            final int[] mag01 = this.mag01; // locals are slightly faster 
	    
	    for (kk = 0; kk < N - M; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+M] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    for (; kk < N-1; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+(M-N)] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
	    mt[N-1] = mt[M-1] ^ (y >>> 1) ^ mag01[y & 0x1];

	    mti = 0;
	    }
  
	y = mt[mti++];
	y ^= y >>> 11;                          // TEMPERING_SHIFT_U(y)
	y ^= (y << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(y)
	y ^= (y << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(y)
	y ^= (y >>> 18);                        // TEMPERING_SHIFT_L(y)

	if (mti >= N)   // generate N words at one time
	    {
	    int kk;
            final int[] mt = this.mt; // locals are slightly faster 
            final int[] mag01 = this.mag01; // locals are slightly faster 
	    
	    for (kk = 0; kk < N - M; kk++)
		{
		z = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+M] ^ (z >>> 1) ^ mag01[z & 0x1];
		}
	    for (; kk < N-1; kk++)
		{
		z = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+(M-N)] ^ (z >>> 1) ^ mag01[z & 0x1];
		}
	    z = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
	    mt[N-1] = mt[M-1] ^ (z >>> 1) ^ mag01[z & 0x1];
	    
	    mti = 0;
	    }
	
	z = mt[mti++];
	z ^= z >>> 11;                          // TEMPERING_SHIFT_U(z)
	z ^= (z << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(z)
	z ^= (z << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(z)
	z ^= (z >>> 18);                        // TEMPERING_SHIFT_L(z)
	
	return ((((long)(y >>> 6)) << 27) + (z >>> 5)) / (double)(1L << 53) < probability;
	}


    public final byte nextByte()
	{
	int y;
	
	if (mti >= N)   // generate N words at one time
	    {
	    int kk;
            final int[] mt = this.mt; // locals are slightly faster 
            final int[] mag01 = this.mag01; // locals are slightly faster 
	    
	    for (kk = 0; kk < N - M; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+M] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    for (; kk < N-1; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+(M-N)] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
	    mt[N-1] = mt[M-1] ^ (y >>> 1) ^ mag01[y & 0x1];

	    mti = 0;
	    }
  
	y = mt[mti++];
	y ^= y >>> 11;                          // TEMPERING_SHIFT_U(y)
	y ^= (y << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(y)
	y ^= (y << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(y)
	y ^= (y >>> 18);                        // TEMPERING_SHIFT_L(y)

	return (byte)(y >>> 24);
	}


    public final void nextBytes(byte[] bytes)
	{
	int y;
	
	for (int x=0;x<bytes.length;x++)
	    {
	    if (mti >= N)   // generate N words at one time
		{
		int kk;
                final int[] mt = this.mt; // locals are slightly faster 
                final int[] mag01 = this.mag01; // locals are slightly faster 
		
		for (kk = 0; kk < N - M; kk++)
		    {
		    y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		    mt[kk] = mt[kk+M] ^ (y >>> 1) ^ mag01[y & 0x1];
		    }
		for (; kk < N-1; kk++)
		    {
		    y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		    mt[kk] = mt[kk+(M-N)] ^ (y >>> 1) ^ mag01[y & 0x1];
		    }
		y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[N-1] = mt[M-1] ^ (y >>> 1) ^ mag01[y & 0x1];
		
		mti = 0;
		}
	    
	    y = mt[mti++];
	    y ^= y >>> 11;                          // TEMPERING_SHIFT_U(y)
	    y ^= (y << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(y)
	    y ^= (y << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(y)
	    y ^= (y >>> 18);                        // TEMPERING_SHIFT_L(y)

	    bytes[x] = (byte)(y >>> 24);
	    }
	}


    public final long nextLong()
	{
	int y;
	int z;

	if (mti >= N)   // generate N words at one time
	    {
	    int kk;
            final int[] mt = this.mt; // locals are slightly faster 
            final int[] mag01 = this.mag01; // locals are slightly faster 
	    
	    for (kk = 0; kk < N - M; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+M] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    for (; kk < N-1; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+(M-N)] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
	    mt[N-1] = mt[M-1] ^ (y >>> 1) ^ mag01[y & 0x1];

	    mti = 0;
	    }
  
	y = mt[mti++];
	y ^= y >>> 11;                          // TEMPERING_SHIFT_U(y)
	y ^= (y << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(y)
	y ^= (y << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(y)
	y ^= (y >>> 18);                        // TEMPERING_SHIFT_L(y)

	if (mti >= N)   // generate N words at one time
	    {
	    int kk;
            final int[] mt = this.mt; // locals are slightly faster 
            final int[] mag01 = this.mag01; // locals are slightly faster 
	    
	    for (kk = 0; kk < N - M; kk++)
		{
		z = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+M] ^ (z >>> 1) ^ mag01[z & 0x1];
		}
	    for (; kk < N-1; kk++)
		{
		z = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+(M-N)] ^ (z >>> 1) ^ mag01[z & 0x1];
		}
	    z = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
	    mt[N-1] = mt[M-1] ^ (z >>> 1) ^ mag01[z & 0x1];
	    
	    mti = 0;
	    }
	
	z = mt[mti++];
	z ^= z >>> 11;                          // TEMPERING_SHIFT_U(z)
	z ^= (z << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(z)
	z ^= (z << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(z)
	z ^= (z >>> 18);                        // TEMPERING_SHIFT_L(z)
	
	return (((long)y) << 32) + (long)z;
	}



    /** Returns a long drawn uniformly from 0 to n-1.  Suffice it to say,
	n must be > 0, or an IllegalArgumentException is raised. */
    public final long nextLong(final long n)
	{
	if (n<=0)
	    throw new IllegalArgumentException("n must be positive");
	
	long bits, val;
	do 
	    {
            int y;
            int z;
    
            if (mti >= N)   // generate N words at one time
                {
                int kk;
                final int[] mt = this.mt; // locals are slightly faster 
                final int[] mag01 = this.mag01; // locals are slightly faster 
	    
                for (kk = 0; kk < N - M; kk++)
                    {
                    y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
                    mt[kk] = mt[kk+M] ^ (y >>> 1) ^ mag01[y & 0x1];
                    }
                for (; kk < N-1; kk++)
                    {
                    y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
                    mt[kk] = mt[kk+(M-N)] ^ (y >>> 1) ^ mag01[y & 0x1];
                    }
                y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
                mt[N-1] = mt[M-1] ^ (y >>> 1) ^ mag01[y & 0x1];
    
                mti = 0;
                }
    
            y = mt[mti++];
            y ^= y >>> 11;                          // TEMPERING_SHIFT_U(y)
            y ^= (y << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(y)
            y ^= (y << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(y)
            y ^= (y >>> 18);                        // TEMPERING_SHIFT_L(y)
    
            if (mti >= N)   // generate N words at one time
                {
                int kk;
                final int[] mt = this.mt; // locals are slightly faster 
                final int[] mag01 = this.mag01; // locals are slightly faster 
                
                for (kk = 0; kk < N - M; kk++)
                    {
                    z = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
                    mt[kk] = mt[kk+M] ^ (z >>> 1) ^ mag01[z & 0x1];
                    }
                for (; kk < N-1; kk++)
                    {
                    z = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
                    mt[kk] = mt[kk+(M-N)] ^ (z >>> 1) ^ mag01[z & 0x1];
                    }
                z = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
                mt[N-1] = mt[M-1] ^ (z >>> 1) ^ mag01[z & 0x1];
                
                mti = 0;
                }
            
            z = mt[mti++];
            z ^= z >>> 11;                          // TEMPERING_SHIFT_U(z)
            z ^= (z << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(z)
            z ^= (z << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(z)
            z ^= (z >>> 18);                        // TEMPERING_SHIFT_L(z)
            
            bits = (((((long)y) << 32) + (long)z) >>> 1);
            val = bits % n;
            } while (bits - val + (n-1) < 0);
        return val;
	}

    /** Returns a random double.  1.0 and 0.0 are both valid results. */
    public final double nextDouble()
	{
	int y;
	int z;

	if (mti >= N)   // generate N words at one time
	    {
	    int kk;
            final int[] mt = this.mt; // locals are slightly faster 
            final int[] mag01 = this.mag01; // locals are slightly faster 
	    
	    for (kk = 0; kk < N - M; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+M] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    for (; kk < N-1; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+(M-N)] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
	    mt[N-1] = mt[M-1] ^ (y >>> 1) ^ mag01[y & 0x1];

	    mti = 0;
	    }
  
	y = mt[mti++];
	y ^= y >>> 11;                          // TEMPERING_SHIFT_U(y)
	y ^= (y << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(y)
	y ^= (y << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(y)
	y ^= (y >>> 18);                        // TEMPERING_SHIFT_L(y)

	if (mti >= N)   // generate N words at one time
	    {
	    int kk;
            final int[] mt = this.mt; // locals are slightly faster 
            final int[] mag01 = this.mag01; // locals are slightly faster 
	    
	    for (kk = 0; kk < N - M; kk++)
		{
		z = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+M] ^ (z >>> 1) ^ mag01[z & 0x1];
		}
	    for (; kk < N-1; kk++)
		{
		z = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+(M-N)] ^ (z >>> 1) ^ mag01[z & 0x1];
		}
	    z = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
	    mt[N-1] = mt[M-1] ^ (z >>> 1) ^ mag01[z & 0x1];
	    
	    mti = 0;
	    }
	
	z = mt[mti++];
	z ^= z >>> 11;                          // TEMPERING_SHIFT_U(z)
	z ^= (z << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(z)
	z ^= (z << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(z)
	z ^= (z >>> 18);                        // TEMPERING_SHIFT_L(z)
	
	/* derived from nextDouble documentation in jdk 1.2 docs, see top */
	return ((((long)(y >>> 6)) << 27) + (z >>> 5)) / (double)(1L << 53);
	}





    public final double nextGaussian()
	{
	if (__haveNextNextGaussian)
	    {
	    __haveNextNextGaussian = false;
	    return __nextNextGaussian;
	    } 
	else 
	    {
	    double v1, v2, s;
	    do 
		{ 
		int y;
		int z;
		int a;
		int b;
		    
                if (mti >= N)   // generate N words at one time
                    {
                    int kk;
                    final int[] mt = this.mt; // locals are slightly faster 
                    final int[] mag01 = this.mag01; // locals are slightly faster 
                    
                    for (kk = 0; kk < N - M; kk++)
                        {
                        y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
                        mt[kk] = mt[kk+M] ^ (y >>> 1) ^ mag01[y & 0x1];
                        }
                    for (; kk < N-1; kk++)
                        {
                        y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
                        mt[kk] = mt[kk+(M-N)] ^ (y >>> 1) ^ mag01[y & 0x1];
                        }
                    y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
                    mt[N-1] = mt[M-1] ^ (y >>> 1) ^ mag01[y & 0x1];
                    
                    mti = 0;
                    }
		
		y = mt[mti++];
		y ^= y >>> 11;                          // TEMPERING_SHIFT_U(y)
		y ^= (y << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(y)
		y ^= (y << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(y)
		y ^= (y >>> 18);                        // TEMPERING_SHIFT_L(y)
		
		if (mti >= N)   // generate N words at one time
		    {
		    int kk;
                    final int[] mt = this.mt; // locals are slightly faster 
                    final int[] mag01 = this.mag01; // locals are slightly faster 
		    
		    for (kk = 0; kk < N - M; kk++)
			{
			z = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
			mt[kk] = mt[kk+M] ^ (z >>> 1) ^ mag01[z & 0x1];
			}
		    for (; kk < N-1; kk++)
			{
			z = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
			mt[kk] = mt[kk+(M-N)] ^ (z >>> 1) ^ mag01[z & 0x1];
			}
		    z = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		    mt[N-1] = mt[M-1] ^ (z >>> 1) ^ mag01[z & 0x1];
		    
		    mti = 0;
		    }
		
		z = mt[mti++];
		z ^= z >>> 11;                          // TEMPERING_SHIFT_U(z)
		z ^= (z << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(z)
		z ^= (z << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(z)
		z ^= (z >>> 18);                        // TEMPERING_SHIFT_L(z)
		
		if (mti >= N)   // generate N words at one time
		    {
		    int kk;
                    final int[] mt = this.mt; // locals are slightly faster 
                    final int[] mag01 = this.mag01; // locals are slightly faster 
		    
		    for (kk = 0; kk < N - M; kk++)
			{
			a = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
			mt[kk] = mt[kk+M] ^ (a >>> 1) ^ mag01[a & 0x1];
			}
		    for (; kk < N-1; kk++)
			{
			a = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
			mt[kk] = mt[kk+(M-N)] ^ (a >>> 1) ^ mag01[a & 0x1];
			}
		    a = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		    mt[N-1] = mt[M-1] ^ (a >>> 1) ^ mag01[a & 0x1];
		    
		    mti = 0;
		    }
		
		a = mt[mti++];
		a ^= a >>> 11;                          // TEMPERING_SHIFT_U(a)
		a ^= (a << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(a)
		a ^= (a << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(a)
		a ^= (a >>> 18);                        // TEMPERING_SHIFT_L(a)
		
		if (mti >= N)   // generate N words at one time
		    {
		    int kk;
                    final int[] mt = this.mt; // locals are slightly faster 
                    final int[] mag01 = this.mag01; // locals are slightly faster 
		    
		    for (kk = 0; kk < N - M; kk++)
			{
			b = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
			mt[kk] = mt[kk+M] ^ (b >>> 1) ^ mag01[b & 0x1];
			}
		    for (; kk < N-1; kk++)
			{
			b = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
			mt[kk] = mt[kk+(M-N)] ^ (b >>> 1) ^ mag01[b & 0x1];
			}
		    b = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		    mt[N-1] = mt[M-1] ^ (b >>> 1) ^ mag01[b & 0x1];
		    
		    mti = 0;
		    }
		
		b = mt[mti++];
		b ^= b >>> 11;                          // TEMPERING_SHIFT_U(b)
		b ^= (b << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(b)
		b ^= (b << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(b)
		b ^= (b >>> 18);                        // TEMPERING_SHIFT_L(b)
		
		/* derived from nextDouble documentation in jdk 1.2 docs, see top */
		v1 = 2 *
		    (((((long)(y >>> 6)) << 27) + (z >>> 5)) / (double)(1L << 53))
		    - 1;
		v2 = 2 * (((((long)(a >>> 6)) << 27) + (b >>> 5)) / (double)(1L << 53))
		    - 1;
		s = v1 * v1 + v2 * v2;
		} while (s >= 1 || s==0);
	    double multiplier = /*Strict*/Math.sqrt(-2 * /*Strict*/Math.log(s)/s);
	    __nextNextGaussian = v2 * multiplier;
	    __haveNextNextGaussian = true;
	    return v1 * multiplier;
	    }
	}
    
    
    
    


    public final float nextFloat()
	{
	int y;
	
	if (mti >= N)   // generate N words at one time
	    {
	    int kk;
            final int[] mt = this.mt; // locals are slightly faster 
            final int[] mag01 = this.mag01; // locals are slightly faster 
	    
	    for (kk = 0; kk < N - M; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+M] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    for (; kk < N-1; kk++)
		{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+(M-N)] ^ (y >>> 1) ^ mag01[y & 0x1];
		}
	    y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
	    mt[N-1] = mt[M-1] ^ (y >>> 1) ^ mag01[y & 0x1];

	    mti = 0;
	    }
  
	y = mt[mti++];
	y ^= y >>> 11;                          // TEMPERING_SHIFT_U(y)
	y ^= (y << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(y)
	y ^= (y << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(y)
	y ^= (y >>> 18);                        // TEMPERING_SHIFT_L(y)

	return (y >>> 8) / ((float)(1 << 24));
	}



    /** Returns an integer drawn uniformly from 0 to n-1.   */
    public final int nextInt(final int n)
	{
	if (n<=0)
	    throw new IllegalArgumentException("n must be positive");
	
	if ((n & -n) == n)  // i.e., n is a power of 2
	    {
	    int y;
	
	    if (mti >= N)   // generate N words at one time
		{
		int kk;
                final int[] mt = this.mt; // locals are slightly faster 
                final int[] mag01 = this.mag01; // locals are slightly faster 
		
		for (kk = 0; kk < N - M; kk++)
		    {
		    y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		    mt[kk] = mt[kk+M] ^ (y >>> 1) ^ mag01[y & 0x1];
		    }
		for (; kk < N-1; kk++)
		    {
		    y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		    mt[kk] = mt[kk+(M-N)] ^ (y >>> 1) ^ mag01[y & 0x1];
		    }
		y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[N-1] = mt[M-1] ^ (y >>> 1) ^ mag01[y & 0x1];
		
		mti = 0;
		}
	    
	    y = mt[mti++];
	    y ^= y >>> 11;                          // TEMPERING_SHIFT_U(y)
	    y ^= (y << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(y)
	    y ^= (y << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(y)
	    y ^= (y >>> 18);                        // TEMPERING_SHIFT_L(y)
	    
	    return (int)((n * (long) (y >>> 1) ) >> 31);
	    }
	
	int bits, val;
	do 
	    {
	    int y;
	    
	    if (mti >= N)   // generate N words at one time
		{
		int kk;
                final int[] mt = this.mt; // locals are slightly faster 
                final int[] mag01 = this.mag01; // locals are slightly faster 
		
		for (kk = 0; kk < N - M; kk++)
		    {
		    y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		    mt[kk] = mt[kk+M] ^ (y >>> 1) ^ mag01[y & 0x1];
		    }
		for (; kk < N-1; kk++)
		    {
		    y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		    mt[kk] = mt[kk+(M-N)] ^ (y >>> 1) ^ mag01[y & 0x1];
		    }
		y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[N-1] = mt[M-1] ^ (y >>> 1) ^ mag01[y & 0x1];
		
		mti = 0;
		}
	    
	    y = mt[mti++];
	    y ^= y >>> 11;                          // TEMPERING_SHIFT_U(y)
	    y ^= (y << 7) & TEMPERING_MASK_B;       // TEMPERING_SHIFT_S(y)
	    y ^= (y << 15) & TEMPERING_MASK_C;      // TEMPERING_SHIFT_T(y)
	    y ^= (y >>> 18);                        // TEMPERING_SHIFT_L(y)
	
	    bits = (y >>> 1);
	    val = bits % n;
	    } while(bits - val + (n-1) < 0);
	return val;
	}
    

    /**
     * Tests the code.
     */
    public static void main(String args[])
        { 
	int j;

	MersenneTwisterFast r;

        // CORRECTNESS TEST
        // COMPARE WITH http://www.math.keio.ac.jp/matumoto/CODES/MT2002/mt19937ar.out
	
	r = new MersenneTwisterFast(new int[]{0x123, 0x234, 0x345, 0x456});
	System.out.println("Output of MersenneTwisterFast with new (2002/1/26) seeding mechanism");
	for (j=0;j<1000;j++)
	    {
	    // first, convert the int from signed to "unsigned"
	    long l = (long)r.nextInt();
	    if (l < 0 ) l += 4294967296L;  // max int value
	    String s = String.valueOf(l);
	    while(s.length() < 10) s = " " + s;  // buffer
	    System.out.print(s + " ");
	    if (j%5==4) System.out.println();	    
	    }

	// SPEED TEST

	final long SEED = 4357;

	int xx; long ms;
	System.out.println("\nTime to test grabbing 100000000 ints");
          
	Random rr = new Random(SEED);
	xx = 0;
	ms = System.currentTimeMillis();
	for (j = 0; j < 100000000; j++)
	    xx += rr.nextInt();
	System.out.println("java.util.Random: " + (System.currentTimeMillis()-ms) + "          Ignore this: " + xx);
	
	r = new MersenneTwisterFast(SEED);
	ms = System.currentTimeMillis();
	xx=0;
	for (j = 0; j < 100000000; j++)
	    xx += r.nextInt();
	System.out.println("Mersenne Twister Fast: " + (System.currentTimeMillis()-ms) + "          Ignore this: " + xx);
	
	// TEST TO COMPARE TYPE CONVERSION BETWEEN
	// MersenneTwisterFast.java AND MersenneTwister.java
          
	System.out.println("\nGrab the first 1000 booleans");
	r = new MersenneTwisterFast(SEED);
	for (j = 0; j < 1000; j++)
	    {
	    System.out.print(r.nextBoolean() + " ");
	    if (j%8==7) System.out.println();
	    }
	if (!(j%8==7)) System.out.println();
	  
	System.out.println("\nGrab 1000 booleans of increasing probability using nextBoolean(double)");
	r = new MersenneTwisterFast(SEED);
	for (j = 0; j < 1000; j++)
	    {
	    System.out.print(r.nextBoolean((double)(j/999.0)) + " ");
	    if (j%8==7) System.out.println();
	    }
	if (!(j%8==7)) System.out.println();
	  
	System.out.println("\nGrab 1000 booleans of increasing probability using nextBoolean(float)");
	r = new MersenneTwisterFast(SEED);
	for (j = 0; j < 1000; j++)
	    {
	    System.out.print(r.nextBoolean((float)(j/999.0f)) + " ");
	    if (j%8==7) System.out.println();
	    }
	if (!(j%8==7)) System.out.println();
	  
	byte[] bytes = new byte[1000];
	System.out.println("\nGrab the first 1000 bytes using nextBytes");
	r = new MersenneTwisterFast(SEED);
	r.nextBytes(bytes);
	for (j = 0; j < 1000; j++)
	    {
	    System.out.print(bytes[j] + " ");
	    if (j%16==15) System.out.println();
	    }
	if (!(j%16==15)) System.out.println();
	
	byte b;
	System.out.println("\nGrab the first 1000 bytes -- must be same as nextBytes");
	r = new MersenneTwisterFast(SEED);
	for (j = 0; j < 1000; j++)
	    {
	    System.out.print((b = r.nextByte()) + " ");
	    if (b!=bytes[j]) System.out.print("BAD ");
	    if (j%16==15) System.out.println();
	    }
	if (!(j%16==15)) System.out.println();

	System.out.println("\nGrab the first 1000 shorts");
	r = new MersenneTwisterFast(SEED);
	for (j = 0; j < 1000; j++)
	    {
	    System.out.print(r.nextShort() + " ");
	    if (j%8==7) System.out.println();
	    }
	if (!(j%8==7)) System.out.println();

	System.out.println("\nGrab the first 1000 ints");
	r = new MersenneTwisterFast(SEED);
	for (j = 0; j < 1000; j++)
	    {
	    System.out.print(r.nextInt() + " ");
	    if (j%4==3) System.out.println();
	    }
	if (!(j%4==3)) System.out.println();

	System.out.println("\nGrab the first 1000 ints of different sizes");
	r = new MersenneTwisterFast(SEED);
	int max = 1;
	for (j = 0; j < 1000; j++)
	    {
	    System.out.print(r.nextInt(max) + " ");
	    max *= 2;
	    if (max <= 0) max = 1;
	    if (j%4==3) System.out.println();
	    }
	if (!(j%4==3)) System.out.println();

	System.out.println("\nGrab the first 1000 longs");
	r = new MersenneTwisterFast(SEED);
	for (j = 0; j < 1000; j++)
	    {
	    System.out.print(r.nextLong() + " ");
	    if (j%3==2) System.out.println();
	    }
	if (!(j%3==2)) System.out.println();

	System.out.println("\nGrab the first 1000 longs of different sizes");
	r = new MersenneTwisterFast(SEED);
	long max2 = 1;
	for (j = 0; j < 1000; j++)
	    {
	    System.out.print(r.nextLong(max2) + " ");
	    max2 *= 2;
	    if (max2 <= 0) max2 = 1;
	    if (j%4==3) System.out.println();
	    }
	if (!(j%4==3)) System.out.println();
          
	System.out.println("\nGrab the first 1000 floats");
	r = new MersenneTwisterFast(SEED);
	for (j = 0; j < 1000; j++)
	    {
	    System.out.print(r.nextFloat() + " ");
	    if (j%4==3) System.out.println();
	    }
	if (!(j%4==3)) System.out.println();

	System.out.println("\nGrab the first 1000 doubles");
	r = new MersenneTwisterFast(SEED);
	for (j = 0; j < 1000; j++)
	    {
	    System.out.print(r.nextDouble() + " ");
	    if (j%3==2) System.out.println();
	    }
	if (!(j%3==2)) System.out.println();

	System.out.println("\nGrab the first 1000 gaussian doubles");
	r = new MersenneTwisterFast(SEED);
	for (j = 0; j < 1000; j++)
	    {
	    System.out.print(r.nextGaussian() + " ");
	    if (j%3==2) System.out.println();
	    }
	if (!(j%3==2)) System.out.println();
	
        }
    }
