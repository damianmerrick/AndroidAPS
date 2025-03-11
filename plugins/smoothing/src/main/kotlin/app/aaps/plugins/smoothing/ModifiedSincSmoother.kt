import kotlin.math.ceil
import kotlin.math.cos
import kotlin.math.exp
import kotlin.math.min
import kotlin.math.pow
import kotlin.math.sin

/**
 * Note: Converted from Java to Kotlin in Android Studio, 3 Mar 2024.
 *
 * A Java implementation of smoothing by a modified sinc
 * kernel (MS or MS1), as described in M. Schmid and U. Diebold,
 * 'Why and how Savitzky-Golay filters should be replaced'
 * The term 'degree' is defined in analogy to Savitzky-Golay (SG) filters;
 * the current MS filters have a similar frequency response as SG filters
 * of the same degree (2, 4, ... 10).
 * <h3>
 * Copyright notice
</h3> *
 * This code is licensed under GNU General Public License (GPLv3) and
 * Creative Commons Attribution-ShareAlike 4.0 (CC-BY-SA 4.0).
 * When using and/or modifying this program for scientific work
 * and the paper on it has been published, please cite the paper:
 * M. Schmid, D. Rath and U. Diebold,
 * 'Why and how Savitzky-Golay filters should be replaced',
 * ACS Measurement Science Au, 2022
 *
 *
 * Author: Michael Schmid, IAP/TU Wien, 2021.
 * https://www.iap.tuwien.ac.at/www/surface/group/schmid
 */
class ModifiedSincSmoother
/**
 * Creates a ModifiedSincSmoother with given degree and given kernel size.
 * This constructor is useful for repeated smoothing operations with
 * the same parameter; then the non-statioc `smooth(double[], double[])`
 * can be used without calculating the kernel and boundary fit weights each
 * time.
 * Otherwise the static `smooth(double[], int, int)`
 * method is more convenient.
 * @param isMS1 if true, uses the MS1 variant, which has a smaller kernel size,
 * at the cost of reduced stopband suppression and more gradual cutoff
 * for degree=2. Otherwise, standard MS kernels are used.
 * @param degree Degree of the filter, must be 2, 4, ... MAX_DEGREE.
 * As for Savitzky-Golay filters, higher degree results in
 * a sharper cutoff in the frequency domain.
 * @param m The half-width of the kernel, must be larger than degree/2.
 * The kernel size is `2*m + 1`.
 */(
    /** Whether MS1 (not MS) filtering should be used  */
    private var isMS1: Boolean,
    /** The degree (2, 4, ... 10)  */
    private var degree: Int, m: Int
) {

    /** The kernel for filtering  */
    private var kernel: DoubleArray

    /** The weights for linear fitting for extending the data at the boundaries  */
    private var fitWeights: DoubleArray

    init {
        require(!(degree < 2 || degree > MAX_DEGREE || degree and 0x1 != 0)) { "Invalid degree " + degree + "; only 2, 4, ... " + MAX_DEGREE + " supported" }
        val mMin = if (isMS1) degree / 2 + 1 else degree / 2 + 2
        require(!(m < mMin)) { "Invalid kernel half-width $m; must be >= $mMin" }
        kernel = makeKernel(isMS1, degree, m)
        fitWeights = makeFitWeights(isMS1, degree, m)
    }

    /**
     * Smooths the data with the parameters passed with the constructor,
     * except for the near-end points.
     * @param data The input data.
     * @param output  The output array; may be null. If `out` is
     * supplied, has the correct size, and is not the input array,
     * it is used for the output.
     * @return  The smoothed data. If `out` is non-null and has
     * the correct size, this is the `out` array.
     * Values within `m` points from boundaries, where
     * the convolution is undefined remain 0 (or retain the
     * previous value, if the supplied `out` array
     * is used).
     */
    fun smoothExceptBoundaries(data: DoubleArray, output: DoubleArray?): DoubleArray {
        var out = output
        if (out == null || out.size != data.size) out = DoubleArray(data.size)
        val radius = kernel.size - 1 //how many additional points we need
        for (i in radius until data.size - radius) {
            var sum = kernel[0] * data[i]
            for (j in 1 until kernel.size) {
                sum += kernel[j] * (data[i - j] + data[i + j])
            }
            out[i] = sum
        }
        return out
    }

    /**
     * Smooths the data with the parameters passed with the constructor,
     * including the near-boundary points. The near-boundary points are
     * handled by weighted linear extrapolation of the data before smoothing.
     * @param data The input data.
     * @param output  The output array; may be null. If `out` is
     * supplied and has the correct size, it is used for the output.
     * @return  The smoothed data. If `out` is non-null and has
     * the correct size, this is the `out` array.
     */
    fun smooth(data: DoubleArray, output: DoubleArray?): DoubleArray {
        var out = output
        val radius = kernel.size - 1
        val extendedData = extendData(data, radius)
        val extendedSmoothed = smoothExceptBoundaries(extendedData, null)
        if (out == null || out.size != data.size) out = DoubleArray(data.size)
        System.arraycopy(extendedSmoothed, radius, out, 0, data.size)
        return out
    }

    /**
     * Extends the data by a weighted fit to a linear function (linear regression).
     * @param data The input data
     * @param m The halfwidth of the kernel. The number of data points
     * contributing to one output value is `2*m + 1`.
     * @return The input data with extrapolated values appended at both ends.
     * At each end, `m` extrapolated points are appended.
     */
    fun extendData(data: DoubleArray, m: Int): DoubleArray {
        val extendedData = DoubleArray(data.size + 2 * m)
        System.arraycopy(data, 0, extendedData, m, data.size)
        val linreg = LinearRegression()
        // linear fit of first points and extrapolate
        val fitLength = min(fitWeights.size.toDouble(), data.size.toDouble()).toInt()
        for (p in 0 until fitLength) linreg.addPointW(p.toDouble(), data[p], fitWeights[p])
        var offset = linreg.getLROffset()
        var slope = linreg.getLRSlope()
        for (p in -1 downTo -m) extendedData[m + p] = offset + slope * p
        // linear fit of last points and extrapolate
        linreg.clear()
        for (p in 0 until fitLength) linreg.addPointW(p.toDouble(), data[data.size - 1 - p], fitWeights[p])
        offset = linreg.getLROffset()
        slope = linreg.getLRSlope()
        for (p in -1 downTo -m) extendedData[(data.size + m) - 1 - p] = offset + slope * p
        return extendedData
    }

    companion object {
        /** This implementation is for a maximum degree of 10  */
        const val MAX_DEGREE = 10

        /** Coefficients for the MS filters, for obtaining a flat passband.
         * The innermost arrays contain a, b, c for the fit
         * kappa = a + b/(c - m)  */
        val CORRECTION_DATA = arrayOf(
            null,  //not defined for degree 0
            null,  //no correction required for degree 2
            null, arrayOf(doubleArrayOf(0.001717576, 0.02437382, 1.64375)), arrayOf(doubleArrayOf(0.0043993373, 0.088211164, 2.359375), doubleArrayOf(0.006146815, 0.024715371, 3.6359375)), arrayOf(
                doubleArrayOf(0.0011840032, 0.04219344, 2.746875), doubleArrayOf(0.0036718843, 0.12780383, 2.7703125)
            )
        )

        /** Coefficients for the MS1 filters, for obtaining a flat passband.
         * The innermost arrays contain a, b, c for the fit
         * kappa = a + b/(c - m)  */
        val CORRECTION_DATA1 = arrayOf(
            null,  //not defined for degree 0
            null,
            arrayOf(doubleArrayOf(0.021944195, 0.050284006, 0.765625)),
            arrayOf(doubleArrayOf(0.0018977303, 0.008476806, 1.2625), doubleArrayOf(0.023064667, 0.13047926, 1.2265625)),
            arrayOf(
                doubleArrayOf(0.0065903002, 0.057929456, 1.915625), doubleArrayOf(0.0023234477, 0.010298849, 2.2726562), doubleArrayOf(0.021046653, 0.16646601, 1.98125)
            ),
            arrayOf(
                doubleArrayOf(9.749618E-4, 0.0020742896, 3.74375),
                doubleArrayOf(0.008975366, 0.09902466, 2.7078125),
                doubleArrayOf(0.0024195414, 0.010064855, 3.296875),
                doubleArrayOf(0.019185117, 0.18953617, 2.784961)
            )
        )

        /**
         * Smooths the data and with the given parameters.
         * When smoothing multiple data sets with the same parameters,
         * using the constructor and then smooth(double[], double[])
         * will be more efficient.
         * @param data   The input data.
         * @param isMS1 if true, uses the MS1 variant, which has a smaller kernel size,
         * at the cost of reduced stopband suppression and more gradual cutoff
         * for degree=2. Otherwise, standard MS kernels are used.
         * @param degree Degree of the filter, must be 2, 4, 6, ... MAX_DEGREE.
         * As for Savitzky-Golay filters, higher degree results in
         * a sharper cutoff in the frequency domain.
         * @param m The half-width of the kernel. The kernel size is
         * `2*m + 1`.
         * @return The smoothed data. Values within `m` points from
         * boundaries, where the convolution is undefined, are set to 0.
         */
        fun smooth(data: DoubleArray, isMS1: Boolean, degree: Int, m: Int): DoubleArray {
            val smoother = ModifiedSincSmoother(isMS1, degree, m)
            return smoother.smooth(data, null)
        }

        /**
         * Calculates the kernel halfwidth m best suited for obtaining a given noise gain.
         * @param isMS1 if true, calculates for the MS1 variant, which has a smaller kernel size,
         * at the cost of reduced stopband suppression and more gradual cutoff
         * for degree=2. Otherwise, standard MS kernels are used.
         * @param degree    The degree n of the kernel
         * @param noiseGain The factor by which white noise should be suppressed.
         * @return The kernel halfwidth m required.
         */
        fun noiseGainToM(isMS1: Boolean, degree: Int, noiseGain: Double): Int {
            val invNoiseGainSqr = 1.0 / (noiseGain * noiseGain)
            val exponent = -2.5 - 0.8 * degree
            val m: Double =
                if (isMS1) -1 + invNoiseGainSqr * (0.543 + 0.4974 * degree) + 0.47 * invNoiseGainSqr.pow(exponent) else -1 + invNoiseGainSqr * (1.494 + 0.4965 * degree) + 0.52 *
                    invNoiseGainSqr
                    .pow(
                    exponent
                )
            return Math.round(m).toInt()
        }

        /**
         * Creates a kernel and returns it.
         * @param isMS1  if true, calculates the kernel for the MS1 variant.
         * Otherwise, standard MS kernels are used.
         * @param degree The degree n of the kernel
         * @param m The half-width of the SG kernel. The kernel size of the
         * filter is `2*m + 1`.
         */
        fun makeKernel(isMS1: Boolean, degree: Int, m: Int): DoubleArray {
            val coeffs = getCoefficients(isMS1, degree, m)
            return makeKernel(isMS1, degree, m, coeffs)
        }

        /**
         * Creates a kernel and returns it.
         * @param isMS1  If true, calculates the kernel for the MS1 variant.
         * Otherwise, standard MS kernels are used.
         * @param degree The degree n of the kernel, i.e. the polynomial degree of a Savitzky-Golay filter
         * with similar passband, must be 2, 4, ... MAX_DEGREE.
         * @param m The half-width of the kernel (the resulting kernel
         * has `2*m+1` elements).
         * @param coeffs Correction parameters for a flatter passband, or
         * null for no correction (used for degree 2).
         * @return  One side of the kernel, starting with the element at the
         * center. Since the kernel is symmetric, only one side with
         * `m+1` elements is needed.
         */
        fun makeKernel(isMS1: Boolean, degree: Int, m: Int, coeffs: DoubleArray?): DoubleArray {
            if (degree < 2 || degree > MAX_DEGREE || degree and 0x01 != 0) throw IllegalArgumentException("Unsupported degree $degree")
            val kernel = DoubleArray(m + 1)
            val nCoeffs = coeffs?.size ?: 0
            var sum = 0.0
            for (i in 0..m) {
                val x = i * (1.0 / (m + 1)) //0 at center, 1 at zero
                val sincArg = Math.PI * 0.5 * (if (isMS1) degree + 2 else degree + 4) * x
                var k: Double = if (i == 0) 1.0 else sin(sincArg) / sincArg
                for (j in 0 until nCoeffs) {
                    k += if (isMS1) coeffs!![j] * x * sin((j + 1) * Math.PI * x) //shorter kernel version, needs more correction terms
                    else {
                        val nu = if (degree / 2 and 0x1 == 0) 2 else 1 //start at 1 for degree 6, 10; at 2 for degree 8
                        coeffs!![j] * x * sin((2 * j + nu) * Math.PI * x)
                    }
                }
                val decay = (if (isMS1) 2 else 4).toDouble() //decay alpha =2: 13.5% at end without correction, 2sqrt2 sigma
                k *= exp(-x * x * decay) + exp(-(x - 2) * (x - 2) * decay) + exp(-(x + 2) * (x + 2) * decay) - 2 * exp(-decay) - exp(-9 * decay)
                kernel[i] = k
                sum += k
                if (i > 0) sum += k //off-center kernel elements appear twice
            }
            for (i in 0..m) kernel[i] *= 1.0 / sum //normalize the kernel to sum = 1
            return kernel
        }

        /**
         * Returns the correction coefficients for a Sinc*Gaussian kernel
         * to flatten the passband.
         * @param isMS1  If true, returns the coefficients for the MS1 variant.
         * Otherwise, coefficients for the standard MS kernels returned.
         * @param degree The polynomial degree of a Savitzky-Golay filter
         * with similar passband, must be 2, 4, ... MAX_DEGREE.
         * @param m The half-width of the kernel.
         * @return  Coefficients z for the x*sin((j+1)*PI*x) terms, or null
         * if no correction is required.
         */
        fun getCoefficients(isMS1: Boolean, degree: Int, m: Int): DoubleArray? {
            val correctionData = if (isMS1) CORRECTION_DATA1 else CORRECTION_DATA
            val corrForDeg = correctionData.get(degree / 2) ?: return null
            val coeffs = DoubleArray(corrForDeg.size)
            for (i in corrForDeg.indices) {
                val abc = corrForDeg[i] // a...c of equation a + b/(c - m)^3
                val cm = abc[2] - m // c - m
                coeffs[i] = abc[0] + abc[1] / (cm * cm * cm)
            }
            return coeffs
        }

        /**
         * Returns the weights for the linear fit used for linear extrapolation
         * at the end. The weight function is a Hann (cos^2) function. For beta=1
         * (the beta value for n=4), it decays to zero at the position of the
         * first zero of the sinc function in the kernel. Larger beta values lead
         * to stronger noise suppression near the edges, but the smoothed curve
         * does not follow the input as well as for lower beta (for high degrees,
         * also leading to more ringing near the boundaries).
         * @param isMS1  if true, returns weights for the MS1 variant.
         * Otherwise, returns weights for the standard MS kernels.
         * @param degree The polynomial degree of a Savitzky-Golay filter
         * with similar passband, must be 2, 4, ... MAX_DEGREE.
         * @param m The half-width of the kernel (the resulting kernel
         * has 2*m+1 elements).
         * @return  The fit weights, with array element [0] corresponding
         * to the data value at the very end.
         */
        fun makeFitWeights(isMS1: Boolean, degree: Int, m: Int): DoubleArray {
            val firstZero = if (isMS1) //the first zero of the sinc function
                (m + 1) / (1 + 0.5 * degree) else (m + 1) / (1.5 + 0.5 * degree)
            val beta = if (isMS1) 0.65 + 0.35 * exp(-0.55 * (degree - 4)) else 0.70 + 0.14 * exp(-0.60 * (degree - 4))
            val fitLength = ceil(firstZero * beta).toInt()
            val weights = DoubleArray(fitLength)
            for (p in 0 until fitLength) weights[p] = sqr(cos(0.5 * Math.PI / (firstZero * beta) * p))
            return weights
        }

        /** Returns the square of a number  */
        fun sqr(x: Double): Double {
            return x * x
        }
    }
}

/** Linear regression is used for extrapolating the data at the boundaries  */
internal class LinearRegression {

    /** sum of weights (number of points if all weights are one  */
    protected var sumWeights = 0.0

    /** sum of all x values  */
    protected var sumX = 0.0

    /** sum of all y values  */
    protected var sumY = 0.0

    /** sum of all x*y products  */
    protected var sumXY = 0.0

    /** sum of all squares of x  */
    protected var sumX2 = 0.0

    /** sum of all squares of y  */
    protected var sumY2 = 0.0

    /** result of the regression: offset  */
    protected var offset = Double.NaN

    /** result of the regression: slope  */
    protected var slope = Double.NaN

    /** whether the results (offset, slope) have been calculated already  */
    protected var calculated = false

    /** Clear the linefit  */
    fun clear() {
        sumWeights = 0.0
        sumX = 0.0
        sumY = 0.0
        sumXY = 0.0
        sumX2 = 0.0
        sumY2 = 0.0
        calculated = false
    }

    /** Add a point x,y with weight  */
    fun addPointW(x: Double, y: Double, weight: Double) {
        sumWeights += weight
        sumX += weight * x
        sumY += weight * y
        sumXY += weight * x * y
        sumX2 += weight * x * x
        sumY2 += weight * y * y
        calculated = false
    }

    /** Do the actual regression calculation  */
    fun calculate() {
        val stdX2TimesN = sumX2 - sumX * sumX * (1 / sumWeights)
        //val stdY2TimesN = sumY2 - sumY * sumY * (1 / sumWeights)
        if (sumWeights > 0) {
            slope = (sumXY - sumX * sumY * (1 / sumWeights)) / stdX2TimesN
            if (java.lang.Double.isNaN(slope)) slope = 0.0 //slope 0 if only one x value
        } else slope = Double.NaN
        offset = (sumY - slope * sumX) / sumWeights
        calculated = true
    }

    /** Returns the offset (intersection on y axis) of the fit  */
    fun getLROffset(): Double {
        if (!calculated) calculate()
        return offset
    }

    /** Returns the slope of the line of the fit  */
    fun getLRSlope(): Double {
        if (!calculated) calculate()
        return slope
    }
}
