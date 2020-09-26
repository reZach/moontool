# moontool
Get information about the moon, in javascript. 

The source of algorithms comes from [Moontool for Windows](http://www.fourmilab.ch/moontoolw/), which the source of the application has been copied to this repo in `/original` for safekeeping/archival.

## Methods available
|Method|Returns|
|---|---|
|getMoonInfo(date)|Returns information about the moon on a given date|
|getEaster(year)|Returns the date of Easter for a given year|


### Examples
```javascript
getMoonInfo(new Date());

// Result
{
    angularDiameterInDeg: 0.5227003811621781,
    distanceInKm: 381017.8168555941,
    distanceToSun: 149979823.99030265​,
    moonAgeInDays: 8.807861270471863​, // min 0, max 29.53
    moonIllumination: 0.6493074320692747​, // 0 new moon, 1 full moon
    moonPhase: 0.298262298998364​, // percentage of the current phase of the moon (new moon (0) > full moon > new moon (1))
    sunAngularDiameter: 0.5317725210369415
}
```

```javascript
getEaster(2020); // -> Sun Apr 12 2020
```

## Usage
This package works both in CJS / ESM formats.
```javascript
const moontool = require("@mrspade/moontool");
moontool.getEaster(2020);

// or

import { getMoonInfo } from "@mrspade/moontool";
getMoonInfo(new Date());
```

## Building locally
`npm i`

## Package
_Ensure you have ran `npm i`_

Run `npm run build`
