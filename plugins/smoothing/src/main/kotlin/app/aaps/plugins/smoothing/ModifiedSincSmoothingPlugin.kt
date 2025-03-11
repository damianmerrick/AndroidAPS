package app.aaps.plugins.smoothing

import ModifiedSincSmoother
import app.aaps.core.data.iob.InMemoryGlucoseValue
import app.aaps.core.data.model.TrendArrow
import app.aaps.core.data.plugin.PluginType
import app.aaps.core.interfaces.logging.AAPSLogger
import app.aaps.core.interfaces.plugin.PluginBase
import app.aaps.core.interfaces.plugin.PluginDescription
import app.aaps.core.interfaces.resources.ResourceHelper
import app.aaps.core.interfaces.sharedPreferences.SP
import app.aaps.core.interfaces.smoothing.Smoothing
import org.json.JSONException
import org.json.JSONObject
import javax.inject.Inject
import javax.inject.Singleton
import kotlin.math.max
import kotlin.math.round

@Singleton
class ModifiedSincSmoothingPlugin @Inject constructor(
    aapsLogger: AAPSLogger,
    rh: ResourceHelper,
    val sp: SP
) : PluginBase(
    PluginDescription()
        .mainType(PluginType.SMOOTHING)
        .pluginIcon(app.aaps.core.ui.R.drawable.ic_timeline_24)
        .pluginName(R.string.ms_smoothing_name)
        .shortName(R.string.smoothing_shortname)
        .preferencesId(R.xml.pref_smoothing_ms)
        .description(R.string.description_ms_smoothing),
    aapsLogger, rh
), Smoothing {

    override fun smooth(data: MutableList<InMemoryGlucoseValue>): MutableList<InMemoryGlucoseValue> {
        var windowSize = data.size
        val noiseGain = sp.getDouble(R.string.key_noise_gain, 0.4)
        var degree = sp.getInt(R.string.key_filter_degree, 10)

        // Ensure we have a valid degree (2, 4, 6, 8 or 10)
        degree = if(degree <= 2) {
            2
        } else if(degree <= 4) {
            4
        } else if(degree <= 6) {
            6
        } else if(degree <= 8) {
            8
        } else {
            10
        }

        val mMin = degree / 2 + 2
        val mMax = 100000
        var m = ModifiedSincSmoother.noiseGainToM(false, degree, noiseGain)
        if(m < mMin) {
            m = mMin
        }
        if(m > mMax) {
            m = mMax
        }

        // ADJUST SMOOTHING WINDOW TO ONLY INCLUDE VALID READINGS - from exponential smoothing code
        // Valid readings include:
        // - Values that actually exist (windowSize may not be larger than sizeRecords)
        // - Values that come in approx. every 5 min. If the time gap between two readings is larger, this is likely due to a sensor error or warmup of a new sensor.d
        // - Values that are not 38 mg/dl; 38 mg/dl reflects an xDrip error state (according to a comment in determine-basal.js)

        //MP: Adjust smoothing window if database size is smaller than the default value + 1 (+1 because the reading before the oldest reading to be smoothed will be used in the calculations
        if (data.size <= windowSize) { //MP standard smoothing window
            windowSize = (data.size - 1).coerceAtLeast(0) //MP Adjust smoothing window to the size of database if it is smaller than the original window size; -1 to always have at least one older
        // value to compare against as a buffer to prevent app crashes
        }

        //MP: Adjust smoothing window further if a gap in the BG database is detected, e.g. due to sensor errors of sensor swaps, or if 38 mg/dl are reported (xDrip error state)
        for (i in 0 until windowSize) {
            if (round((data[i].timestamp - data[i + 1].timestamp) / (1000.0 * 60)) >= 12) { //MP: 12 min because a missed reading (i.e. readings coming in after 10 min) can occur for various reasons, like walking away from the phone or reinstalling AAPS
                //if (Math.round((data.get(i).date - data.get(i + 1).date) / 60000L) <= 7) { //MP crashes the app, useful for testing
                windowSize =
                    i + 1 //MP: If time difference between two readings exceeds 7 min, adjust windowSize to *include* the more recent reading (i = reading; +1 because windowSize reflects number of valid readings);
                break
            } else if (data[i].value == 38.0) {
                windowSize = i //MP: 38 mg/dl reflects an xDrip error state; Chain of valid readings ends here, *exclude* this value (windowSize = i; i + 1 would include the current value)
                break
            }
        }

        if (windowSize < 4) { //MP: Require a valid windowSize of at least 4 readings
            // Not enough data for smoothing, so just replicate the existing data array
            for (i in 0 until data.size) {
                data[i].smoothed = max(data[i].value, 39.0) // Make 39 the smallest value as smaller values trigger errors (xDrip error state = 38)
                data[i].trendArrow = TrendArrow.NONE
            }

            return data
        }

        // MS smoothing starts here!

        // Construct double array of input values from the valid window of data (reversed so that the oldest value is first)
        val inValues = DoubleArray(windowSize)
        for (i in 0 until windowSize) {
            inValues[i] = data[windowSize - 1 - i].value
        }

        // Smooth the values into another array of output values
        val outValues = ModifiedSincSmoother.smooth(inValues, false, degree,  m)

        // Reverse the output values and use them to populate the smoothed entries of the data array
        for(i in 0 until windowSize) {
            data[i].smoothed = max(outValues[windowSize - 1 - i], 39.0)
            data[i].trendArrow = TrendArrow.NONE
        }

        return data
    }

    override fun configuration(): JSONObject {
        val c = JSONObject()
        try {
            c.put(rh.gs(R.string.key_noise_gain), sp.getDouble(R.string.key_noise_gain, 0.4))
            c.put(rh.gs(R.string.key_filter_degree), sp.getInt(R.string.key_filter_degree, 10))
        } catch (e: JSONException) {
            e.printStackTrace()
        }
        return c
    }

    override fun applyConfiguration(configuration: JSONObject) {
        try {
            if (configuration.has(rh.gs(R.string.key_noise_gain))) sp.putDouble(rh.gs(R.string.key_noise_gain), configuration.getDouble(rh.gs(R.string.key_noise_gain)))
            if (configuration.has(rh.gs(R.string.key_filter_degree))) sp.putInt(rh.gs(R.string.key_filter_degree), configuration.getInt(rh.gs(R.string.key_filter_degree)))
        } catch (e: JSONException) {
            e.printStackTrace()
        }
    }
}