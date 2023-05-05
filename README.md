# flycatcher_alarm_calls_song_discrimination

This repository hosts data and R codes for Bliard L., Qvarnström A., Wheatcroft D. (2021). The role of introductory alarm calls for song discrimination in *Ficedula* flycatchers. Animal Behaviour. https://doi.org/10.1016/j.anbehav.2021.05.018

Data and code also available on Zenodo https://doi.org/10.5281/zenodo.4607842


## GENERAL INFORMATION

1. Title: Data and scripts from "The role of introductory alarm calls for song discrimination in *Ficedula* flycatchers".

2. Author Information:

	A.  Name: Louis Bliard
		Institution: University of Zurich
		Address: Winterthurerstrasse 190, 8057 Zurich, Switzerland
		Email: louis.bliard@uzh.ch / louis.bliard@evobio.eu

        B.  Name: Anna Qvarnstöm
		Institution: Uppsala University
		Address: Norbyvägen 18 D, 75236 Uppsala, Sweden

	C.  Name: David Wheatcroft
		Institution: Stockholm University
		Address: Svante Arrhenius Väg 18B, 11418 Stockholm, Sweden
		Email: david.wheatcroft@zoologi.su.se

3. Date of data collection: 2019-05-03 – 2019-06-21

4. Geographic location of data collection: Öland, Kalmar County, Sweden 

5. Information about funding sources that supported the collection of the data: Vetenskapsrådet grant ‘VR 2016-05138’


## SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: CC-BY 4.0

2. Links to publications that cite or use the data: https://doi.org/10.1016/j.anbehav.2021.05.018

3. Links to other publicly accessible locations of the data: None

4. Links/relationships to ancillary data sets: https://doi.org/10.1038/s41559-017-0192

5. Was data derived from another source? No

6. Recommended citation for this dataset: Bliard Louis, Qvarnström Anna, & Wheatcroft David. (2021). Data and scripts from: The role of introductory alarm calls for song discrimination in Ficedula flycatchers [Data set]. Zenodo. https://doi.org/10.5281/zenodo.4607842


## DATA & FILE OVERVIEW

1. File List: 
- `data_adult_experiment.txt`
- `analysis_adult_experiment.r`
- `data_nestling_experiment.txt`
- `analysis_nestling_experiment.r`
- `data_nestling_wheatcroft2017.txt`
- `analysis_nestling_data_wheatcroft2017.r`
- `comparison_call_or_not.txt`

2. Relationship between files, if important: 

The dataset `data_adult_experiment.txt` was used for the analysis of the experiment on adult males, and analyses can be reproduced using the R script `analysis_adult_experiment.r`.

The dataset `data_nestling_experiment.txt` was used for the analysis of the experiment on nestlings, and analyses can be reproduced using the R script `analysis_nestling_experiment.r`.

The dataset `data_nestling_wheatcroft2017.txt` was used for the re-analysis of the experiment on nestlings from Wheatcroft and Qvarnström data (https://doi.org/10.1038/s41559-017-0192), and analyses can be reproduced using the R script `analysis_nestling_data_wheatcroft2017.r`.

The dataset `comparison_call_or_not.txt` was used for the analysis in Appendix I.
 

## METHODOLOGICAL INFORMATION

1. Description of methods used for collection/generation of data: 

Dataset “data_nestling_experiment.txt”

12-day-old pied flycatcher nestlings were exposed to playbacks of different sound recordings. All the nestlings of a clutch were temporarily removed from their nest and placed on a layer of moss inside an experimental nestbox. Each nestling was marked with symbols on the head using correction fluid, to enable individual recognition of the nestlings in video recordings. The nestbox was set up at least 200 m from the closest known flycatcher nest to avoid hearing adult flycatchers singing or alarm calling. An ‘Eco Extreme Speaker’ (Grace Digital, Poway, CA, U.S.A.) was positioned approximately 1 m from the nestbox and connected to a smartphone. We waited a minimum of 5 min after the nestlings were placed inside the experimental nestbox to start a trial. Each trial proceeded as follows. After 1 min of silence, one of two treatments was played back: collared flycatcher songs, each preceded by a collared flycatcher call, or pied flycatcher songs, each preceded by a collared flycatcher call. Each playback lasted 1 min. Then, a minute of silence was observed again, followed by a recording from the other treatment. Thus, each collared flycatcher nestling was exposed to the two different treatments, alternating the order between each trial. This playback experiment was performed on 222 nestlings from 45 nests. After each trial, nestlings were individually ringed and weighed, and then brought back to their natal nestbox. The whole trial was watched via a live video feed and video recorded using a digital video recorder (PV-1000, Lawmate International, Taipei, Taiwan). Subsequently, the video recording of each experimental trial was analysed by a single observer. We noted the number of times each nestling begged, opened its gape, looked up and shifted its position inside the experimental nestbox.

Dataset “data_adult_experiment.txt”

This experiment was mimicking a territorial intrusion by another male flycatcher. Each responding male was used only once. A trial consisted of a 5 min period of silence where only the speaker was set up and a wooden collared flycatcher dummy was mounted on the nestbox. During the initial silence period, we ensured that males did not respond to the dummy itself to be sure that they were not responding to cues other than the following playback. The dummy was used to make this experiment comparable to earlier studies. This was followed by one of the two experimental treatments for a duration of 10 min. During the 10 min of the playback, we recorded the behavioural response of the focal male, noting for every minute whether the male was alarm calling, whether it was wing flicking and its minimal perching distance to the speaker.

Dataset “data_nestling_wheatcroft2017.txt”

Similar methods as for the “data_nestling_experiment.txt”, see more details in https://doi.org/10.1038/s41559-017-0192 


2. Methods for processing the data: Raw data

3. Instrument- or software-specific information needed to interpret the data: 
R v.3.5.1 https://www.r-project.org/
JAGS (Just Another Gibbs Sampler) https://mcmc-jags.sourceforge.io/

4. People involved with sample collection, processing, analysis and/or submission: Louis Bliard

### DATA-SPECIFIC INFORMATION FOR: `data_nestling_experiment.txt`

1. Number of variables: 21

2. Number of cases/rows: 444

3. Variable List: 
- "Date" = Day when the trial was performed (yyyy-mm-dd)
- "Area" = Patch of forest where the trial was performed
- "Box" = Nestbox where the trial was performed
- "Chick_ID" = Identity of the nestling on the video recording (number of dots on its head)
- "Ring" = Ring number of the nestling
- "Tarsus" = Tarsus length of the nestling (mm)
- "Mass" = Mass of the nestling (g)
- "Order" = Order of the treatment
- "Treatment" = Experimental treatment used for the trial
- "Playback_ID" = Playback file used for the trial
- "Group_size" = Number of nestlings in the nestbox
- "Beg_rest" = Number of time the nestling begged during the minute of silence of the trial
- "Jump_rest" = Number of time the nestling jumped during the minute of silence of the trial
- "Gape_rest" = Number of time the nestling opened its gape during the minute of silence of the trial
- "Look_rest" = Number of time the nestling looked up during the minute of silence of the trial
- "Shift_rest" = Number of time the nestling shifted its position in the nestbox during the minute of silence of the trial
- "Beg_play" = Number of time the nestling begged during the minute of playback of the trial
- "Jump_play" = Number of time the nestling jumped during the minute of playback of the trial
- "Gape_play" = Number of time the nestling opened its gape during the minute of playback of the trial
- "Look_play" = Number of time the nestling looked up during the minute of playback of the trial
- "Shift_play" = Number of time the nestling shifted its position in the nestbox during the minute of playback of the trial

4. Missing data codes: NA

### DATA-SPECIFIC INFORMATION FOR: `data_adult_experiment.txt`

1. Number of variables: 11

2. Number of cases/rows: 210

3. Variable List: 
- "area" = Patch of forest where the trial was performed
- "nestbox" = Nestbox where the trial was performed
- "day" = Day when the trial was performed (yyyy-mm-dd)
- "time" = Time when the trial was performed, ignore the date part (hh:mm:ss)
- "playback_id" = Playback file used for the trial
- "treatment" = Experimental treatment used for the trial
- "minute" = Describe which minute of the trial the line of data refers to. The first 5 minutes of each trial are not reported since the playback was not yet on (see methods)
- "minimum_distance" = Minimal perching distance (in meters) of the bird within this minute of the trial. If the bird was not seen, the minimal distance was set to 40 meters
- "alarm" = Whether the bird alarm called or not during this minute of the trial
- "wing_flick" = Whether the bird was seen wing flicking during this minute of the trial
- "graded_resp" = Behavioural graded response of the bird for this minute of the trial (see methods)

4. Missing data codes: NA

### DATA-SPECIFIC INFORMATION FOR: `data_nestling_wheatcroft2017.txt`

1. Number of variables: 23

2. Number of cases/rows: 235

3. Variable List: 
- "year" = Year the trial was performed
- "area" = Patch of forest where the trial was performed
- "box = Nestbox where the trial was performed
- "nest" = Unique identifier of the nestbox, including the area and the year
- "date" = Day when the trial was performed (m.dd)
- "nestbox" = Nestbox where the trial was performed
- "species" = Species of the nestlings used. It is only collared flycatchers, since only these trials were used in the re-analysis (see methods)
- "chick" = ID of the nestling within its clutch (nestlings are numbered from 1 to X)
- "nestchick" = Unique identifier of the nestling within its nest
- "order" = Order of the treatment
- "playback" = Species of the song played back. It is only collared flycatchers, since only these trials were used in the re-analysis (see methods)
- "file" = File number of the playback sound file
- "recording" = Name of the playback file
- "Songs.with.call" = Number of song phrases starting with a call within the minute of playback
- "Songs.without.call" = Number of song phrases not starting with a call within the minute of playback
- "Call.bin" = Whether at least one of the song phrases started with a call within the minute of playback
- "Call.prop" = Proportion of song phrases starting with a call within the minute of playback
- "source" = Whether the trial was done in an area in sympatry (sym) or allopatry (allo) with pied flcyatchers
- "treatment" = Conspecific, since we only used collared flycatcher nestlings exposed to collared flycatchers songs in the re-analysis (see methods)
- "beg.rest" = Number of time the nestling begged during the minute of silence of the trial
- "beg.play" = Number of time the nestling begged during the minute of playback of the trial
- "clutch.size" = Number of nestlings in the nestbox
- "mass" = Mass of the nestling (g)

4. Missing data codes: NA
