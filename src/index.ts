import { PhaseResult } from "./types";

/* Astronomical constants */
const epoch: number = 2444238.5;                // 1980 January 0.0

/* Constants defining the Sun's apparent orbit */
const elonge: number = 278.833540;              // Ecliptic longitude of the Sun at epoch 1980.0
const elongp: number = 282.596403;              // Ecliptic longitude of the Sun at perigee
const eccent: number = 0.016718;                // Eccentricity of Earth's orbit
const sunsmax: number = 1.495985e8;             // Semi-major axis of Earth's orbit, km
const sunangsiz: number = 0.533128;             // Sun's angular size, degrees, at semi-major axis distance

/* Elements of the Moon's orbit, epoch 1980.0 */
const mmlong: number = 64.975464;               // Moon's mean longitude at the epoch
const mmlongp: number = 349.383063;             // Mean longitude of the perigee at the epoch
const mlnode: number = 151.950429;              // Mean longitude of the node at the epoch
const minc: number = 5.145396;                  // Inclination of the Moon's orbit
const mecc: number = 0.054900;                  // Eccentricity of the Moon's orbit
const mangsiz: number = 0.5181;                 // Moon's angular size at distance a from Earth
const msmax: number = 384401.0;                 // Semi-major axis of Moon's orbit in km
const mparallax: number = 0.9507;               // Parallax at distance a from Earth
const synmonth: number = 29.53058868;           // Synodic month (new Moon to new Moon)
const lunatbase: number = 2423436.0;            // Base date for E. W. Brown's numbered series of lunations (1923 January 16)

/*  Properties of the Earth  */
const earthrad: number = 6378.16;               // Radius of Earth in kilometres
const PI: number = 3.14159265358979323846;      // Assume not near black hole nor in Tennessee

/* More consts */
const epsilon: number = 0.000001;

/*  Handy mathematical functions  */
function sgn(x: number): number { return (((x) < 0) ? -1 : ((x) > 0 ? 1 : 0)) }             // Extract sign
function abs(x: number): number { return ((x) < 0 ? (-(x)) : (x)) }                         // Absolute val
function fixAngle(a: number): number { return ((a) - 360.0 * (Math.floor((a) / 360.0))) }   // Fix angle   
function toRad(d: number): number { return ((d) * (PI / 180.0)) }                           // Deg->Rad    
function toDeg(d: number): number { return ((d) * (180.0 / PI)) }                           // Rad->Deg    
function dsin(x: number): number { return (Math.sin(toRad((x)))) }                          // Sin from deg
function dcos(x: number): number { return (Math.cos(toRad((x)))) }                          // Cos from deg

/**
 * Convert a Date to astronomica Julian time 
 * (i.e. Julian date plus day fraction, expressed as a double).
 * @param date 
 */
function toJulianTime(date: Date): number {

    // Algorithm as given in Meeus, Astronomical Algorithms
    let year: number,
        month: number,
        day: number,
        hour: number,
        minute: number,
        second: number,
        millisecond: number;

    year = date.getFullYear();
    month = date.getMonth() + 1;
    day = date.getDate();
    hour = date.getHours();
    minute = date.getMinutes();
    second = date.getSeconds();
    millisecond = date.getMilliseconds();

    let isjulian = isJulianDate(year, month, day);
    let m: number = month > 2 ? month : month + 12;
    let y: number = month > 2 ? year : year - 1;
    let d: number = day + hour / 24.0 + minute / 1440.0 + (second + millisecond / 1000.0) / 86400.0;
    let b: number = isjulian ? 0 : 2 - y / 100 + y / 100 / 4;

    return Math.floor((365.25 * (y + 4716)) + Math.floor(30.6001 * (m + 1)) + d + b - 1524.5);
}

/**
 * Returns whether a given date is a Julian date or not.
 * https://stackoverflow.com/a/14554483/1837080
 * @param year
 * @param month 
 * @param day 
 */
function isJulianDate(year: number, month: number, day: number): boolean {

    // All dates prior to 1582 are in the Julian calendar
    if (year < 1582)
        return true;
    // All dates after 1582 are in the Gregorian calendar
    else if (year > 1582)
        return false;
    else {
        // If 1582, check before October 4 (Julian) or after October 15 (Gregorian)
        if (month < 10)
            return true;
        else if (month > 10)
            return false;
        else {
            if (day < 5)
                return true;
            else if (day > 14)
                return false;
            else
                throw "Any date in the range 10/5/1582 to 10/14/1582 is invalid!";
        }
    }
}

/**
 * Convert Julian date to year, month, day.
 * @param td 
 * @param yy 
 * @param mm 
 * @param dd 
 */
function jyear(td: number, yy: number, mm: number, dd: number): { year: number, month: number, day: number } {
    let z: number,
        f: number,
        a: number,
        alpha: number,
        b: number,
        c: number,
        d: number,
        e: number;

    td += 0.5;
    z = Math.floor(td);
    f = td - z;

    if (z < 2299161.0) {
        a = z;
    } else {
        alpha = Math.floor((z - 1867216.25) / 36524.25);
        a = z + 1 + alpha - Math.floor(alpha / 4);
    }

    b = a + 1524;
    c = Math.floor((b - 122.1) / 365.25);
    d = Math.floor(365.25 * c);
    e = Math.floor((b - d) / 30.6001);

    return {
        day: Math.floor(b - d - Math.floor(30.6001 * e) + f),
        month: Math.floor((e < 14) ? (e - 1) : (e - 13)),
        year: Math.floor((mm > 2) ? (c - 4716) : (c - 4715))
    };
}

/**
 * Convert Julian time to hour, minutes, and seconds.
 * @param j 
 * @param h 
 * @param m 
 * @param s 
 */
function jhms(j: number): { hour: number, minute: number, second: number } {
    let ij: number;

    j += 0.5;			                                        // Astronomical to civil
    ij = Math.floor(((j - Math.floor(j)) * 86400.0) + 0.5);     // Round to nearest second

    return {
        hour: Math.floor(ij / 3600),
        minute: Math.floor((ij / 60) % 60),
        second: Math.floor(ij % 60)
    };
}

/**
 * Determine day of the week for a given Julian day.
 * @param j 
 */
function jwday(j: number): number {
    return (Math.floor(j + 1.5)) % 7;
}

/**
 * Calculates time of the mean new Moon for a given base date. 
 * This argument K to this function is the precomputed synodic month index, given by:

    K = (year - 1900) * 12.3685

    where year is expressed as a year and fractional year.
 * @param sdate 
 * @param k 
 */
function meanphase(sdate: number, k: number): number {
    let t: number, t2: number, t3: number, nt1: number;

    /* Time in Julian centuries from 1900 January 0.5 */
    t = (sdate - 2415020.0) / 36525;
    t2 = t * t;                       // Square for frequent use
    t3 = t2 * t;                      // Cube for frequent use

    nt1 = 2415020.75933 + synmonth * k
        + 0.0001178 * t2
        - 0.000000155 * t3
        + 0.00033 * dsin(166.56 + 132.87 * t - 0.009173 * t2);

    return nt1;
}

/**
 * Given a K value used to determine the mean phase of the new moon, 
 * and a phase selector (0.0, 0.25, 0.5, 0.75), obtain the true, corrected phase time.
 * @param k 
 * @param phase 
 */
function truephase(k: number, phase: number): number {
    let t: number,
        t2: number,
        t3: number,
        pt: number,
        m: number,
        mprime: number,
        f: number;
    let apcor = false;

    k += phase;                       // Add phase to new moon time
    t = k / 1236.85;                  // Time in Julian centuries from 1900 January 0.5
    t2 = t * t;                       // Square for frequent use
    t3 = t2 * t;                      // Cube for frequent use
    pt = 2415020.75933                // Mean time of phase
        + synmonth * k
        + 0.0001178 * t2
        - 0.000000155 * t3
        + 0.00033 * dsin(166.56 + 132.87 * t - 0.009173 * t2);

    m = 359.2242                      // Sun's mean anomaly
        + 29.10535608 * k
        - 0.0000333 * t2
        - 0.00000347 * t3;
    mprime = 306.0253                 // Moon's mean anomaly
        + 385.81691806 * k
        + 0.0107306 * t2
        + 0.00001236 * t3;
    f = 21.2964                       // Moon's argument of latitude
        + 390.67050646 * k
        - 0.0016528 * t2
        - 0.00000239 * t3;
    if ((phase < 0.01) || (abs(phase - 0.5) < 0.01)) {

        /* Corrections for New and Full Moon */

        pt += (0.1734 - 0.000393 * t) * dsin(m)
            + 0.0021 * dsin(2 * m)
            - 0.4068 * dsin(mprime)
            + 0.0161 * dsin(2 * mprime)
            - 0.0004 * dsin(3 * mprime)
            + 0.0104 * dsin(2 * f)
            - 0.0051 * dsin(m + mprime)
            - 0.0074 * dsin(m - mprime)
            + 0.0004 * dsin(2 * f + m)
            - 0.0004 * dsin(2 * f - m)
            - 0.0006 * dsin(2 * f + mprime)
            + 0.0010 * dsin(2 * f - mprime)
            + 0.0005 * dsin(m + 2 * mprime);
        apcor = true;
    } else if ((abs(phase - 0.25) < 0.01 || (abs(phase - 0.75) < 0.01))) {
        pt += (0.1721 - 0.0004 * t) * dsin(m)
            + 0.0021 * dsin(2 * m)
            - 0.6280 * dsin(mprime)
            + 0.0089 * dsin(2 * mprime)
            - 0.0004 * dsin(3 * mprime)
            + 0.0079 * dsin(2 * f)
            - 0.0119 * dsin(m + mprime)
            - 0.0047 * dsin(m - mprime)
            + 0.0003 * dsin(2 * f + m)
            - 0.0004 * dsin(2 * f - m)
            - 0.0006 * dsin(2 * f + mprime)
            + 0.0021 * dsin(2 * f - mprime)
            + 0.0003 * dsin(m + 2 * mprime)
            + 0.0004 * dsin(m - 2 * mprime)
            - 0.0003 * dsin(2 * m + mprime);
        if (phase < 0.5)
            /* First quarter correction */
            pt += 0.0028 - 0.0004 * dcos(m) + 0.0003 * dcos(mprime);
        else
            /* Last quarter correction */
            pt += -0.0028 + 0.0004 * dcos(m) - 0.0003 * dcos(mprime);
        apcor = true;
    }

    if (!apcor) {
        // exit nicely?
        throw "Error calculating moon phase!";
    }
    return pt;
}

/**
 * Find time of phases of the moon which surround the current date.
 * Five phases are found, starting and ending with the new moons 
 * which bound the current lunation.
 * @param sdate 
 * @param phases 
 */
function phasehunt(sdate: number, phases: number[]): number[] {
    let adate: number, k1: number, k2: number, nt1: number, nt2: number;
    let yy: number, mm: number, dd: number; // THESE ARE INTS!

    adate = sdate - 45;

    let jyearResult = jyear(adate, yy, mm, dd);
    yy = jyearResult.year;
    mm = jyearResult.month;
    dd = jyearResult.day;
    k1 = Math.floor((yy + ((mm - 1) * (1.0 / 12.0)) - 1900) * 12.3685);

    adate = nt1 = meanphase(adate, k1);
    while (true) {
        adate += synmonth;
        k2 = k1 + 1;
        nt2 = meanphase(adate, k2);
        if (nt1 <= sdate && nt2 > sdate)
            break;
        nt1 = nt2;
        k1 = k2;
    }

    phases[0] = truephase(k1, 0.0);
    phases[1] = truephase(k1, 0.25);
    phases[2] = truephase(k1, 0.5);
    phases[3] = truephase(k1, 0.75);
    phases[4] = truephase(k2, 0.0);

    return phases;
}

/** 
 * Solve the equation of Kepler.
 * @param m 
 * @param ecc 
 */
function kepler(m: number, ecc: number): number {
    let e: number, delta: number;

    m = toRad(m);
    e = m;

    do {
        delta = e - ecc * Math.sin(e) - m;
        e -= delta / (1 - ecc * Math.cos(e));
    } while (abs(delta) > epsilon);
    return e;
}

/** 
 * Calculate phase of moon as a fraction.
 * The argument is the time for which the phase is requested,
    expressed as a Julian date and fraction. Returns the terminator
    phase angle as a percentage of a full circle (i.e., 0 to 1), and
    stores into pointer arguments the illuminated fraction of the
    Moon's disc, the Moon's age in days and fraction, the distance of
    the Moon from the centre of the Earth, and the angular diameter
    subtended by the Moon as seen by an observer at the centre of the
    Earth.

 * @param {number} julianDate
 * @return {PhaseResult}
 */
function getMoonPhase(julianDate: number): PhaseResult {
    let Day: number,
        N: number,
        M: number,
        Ec: number,
        Lambdasun: number,
        ml: number,
        MM: number,
        MN: number,
        Ev: number,
        Ae: number,
        A3: number,
        MmP: number,
        mEc: number,
        A4: number,
        lP: number,
        V: number,
        lPP: number,
        NP: number,
        y: number,
        x: number,
        Lambdamoon: number,
        BetaM: number,
        MoonAge: number,
        MoonPhase: number,
        MoonDist: number,
        MoonDFrac: number,
        MoonAng: number,
        MoonPar: number,
        F: number,
        SunDist: number,
        SunAng: number;

    /* Calculation of the Sun's position */
    Day = julianDate - epoch;                                                // Date within epoch
    N = fixAngle((360 / 365.2422) * Day);                               // Mean anomaly of the Sun
    M = fixAngle(N + elonge - elongp);                                  // Convert from perigee co-ordinates to epoch 1980.0

    Ec = kepler(M, eccent);                                             // Solve equation of Kepler
    Ec = Math.sqrt((1 + eccent) / (1 - eccent)) * Math.tan(Ec / 2);
    Ec = 2 * toDeg(Math.atan(Ec));                                      // True anomaly
    Lambdasun = fixAngle(Ec + elongp);                                  // Sun's geocentric ecliptic longitude

    /* Orbital distance factor */
    F = ((1 + eccent * Math.cos(toRad(Ec))) / (1 - eccent * eccent));
    SunDist = sunsmax / F;                                              // Distance to Sun in km
    SunAng = F * sunangsiz;                                             // Sun's angular size in degrees

    /* Calculation of the Moon's position */

    /* Moon's mean longitude */
    ml = fixAngle(13.1763966 * Day + mmlong);

    /* Moon's mean anomaly */
    MM = fixAngle(ml - 0.1114041 * Day - mmlongp);

    /* Moon's ascending node mean longitude */
    MN = fixAngle(mlnode - 0.0529539 * Day);

    /* Evection */
    Ev = 1.2739 * Math.sin(toRad(2 * (ml - Lambdasun) - MM));

    /* Annual equation */
    Ae = 0.1858 * Math.sin(toRad(M));

    /* Correction term */
    A3 = 0.37 * Math.sin(toRad(M));

    /* Corrected anomaly */
    MmP = MM + Ev - Ae - A3;

    /* Correction for the equation of the centre */
    mEc = 6.2886 * Math.sin(toRad(MmP));

    /* Another correction term */
    A4 = 0.214 * Math.sin(toRad(2 * MmP));

    /* Corrected longitude */
    lP = ml + Ev + mEc - Ae + A4;

    /* Variation */
    V = 0.6583 * Math.sin(toRad(2 * (lP - Lambdasun)));

    /* True longitude */
    lPP = lP + V;

    /* Corrected longitude of the node */
    NP = MN - 0.16 * Math.sin(toRad(M));

    /* Y inclination coordinate */
    y = Math.sin(toRad(lPP - NP)) * Math.cos(toRad(minc));

    /* X inclination coordinate */
    x = Math.cos(toRad(lPP - NP));

    /* Ecliptic longitude */
    Lambdamoon = toDeg(Math.atan2(y, x));
    Lambdamoon += NP;

    /* Ecliptic latitude */
    BetaM = toDeg(Math.asin(Math.sin(toRad(lPP - NP)) * Math.sin(toRad(minc))));

    /* Calculation of the phase of the Moon */

    /* Age of the Moon in degrees */
    MoonAge = lPP - Lambdasun;

    /* Phase of the Moon */
    MoonPhase = (1 - Math.cos(toRad(MoonAge))) / 2;

    /* Calculate distance of moon from the centre of the Earth */

    MoonDist = (msmax * (1 - mecc * mecc)) /
        (1 + mecc * Math.cos(toRad(MmP + mEc)));

    /* Calculate Moon's angular diameter */

    MoonDFrac = MoonDist / msmax;
    MoonAng = mangsiz / MoonDFrac;

    /* Calculate Moon's parallax */

    MoonPar = mparallax / MoonDFrac;

    return {
        moonIllumination: MoonPhase,
        moonAgeInDays: synmonth * (fixAngle(MoonAge) / 360.0),
        distanceInKm: MoonDist,
        angularDiameterInDeg: MoonAng,
        distanceToSun: SunDist,
        sunAngularDiameter: SunAng,
        moonPhase: fixAngle(MoonAge) / 360.0
    };
}

/**
 * Get moon information on a given date
 * @param date 
 */
function getMoonInfo(date: Date): PhaseResult {
    if (typeof date === "undefined" || date === null)
        return {
            moonPhase: 0,
            moonIllumination: 0,
            moonAgeInDays: 0,
            distanceInKm: 0,
            angularDiameterInDeg: 0,
            distanceToSun: 0,
            sunAngularDiameter: 0
        };

    return getMoonPhase(toJulianTime(date));
}

/**
 * Return the date of Easter for a given date (ie. year)
 * @param date 
 */
function getEaster(date: Date): Date {
    // Easter is on the first Sunday following a full
    // moon after the vernal equinox

    // Church recognizes the vernal equinox on March 21st;
    // js date object's month is 0-based
    let start = new Date(date.getFullYear(), 2, 21);

    // Continue to calculate moon info, until we find the 
    // first full moon after the vernal equinox
    let fullMoon = start;
    let previousMoonInfo;
    let moonInfo;    
    let gettingDarker = undefined;

    do {
        previousMoonInfo = getMoonInfo(fullMoon);
        fullMoon.setDate(fullMoon.getDate() + 1);
        moonInfo = getMoonInfo(fullMoon);

        // Initially set if we must currently wait for the moon to grow dimmer,
        // before it grows bright again for a full moon
        if (typeof gettingDarker === "undefined"){
            gettingDarker = moonInfo.moonIllumination < previousMoonInfo.moonIllumination;
        } else if (gettingDarker && moonInfo.moonIllumination > previousMoonInfo.moonIllumination){
            
            // Once the moon has finished getting darker,
            // change this variable so we can check that it continues to grow
            // brighter (so we know when we've found our full moon)
            gettingDarker = false;
        }
    } while(gettingDarker && moonInfo.moonIllumination < previousMoonInfo.moonIllumination ||
            !gettingDarker && moonInfo.moonIllumination > previousMoonInfo.moonIllumination);

    // We found a full moon, go back a day since
    // we've gone 1 day too far
    fullMoon.setDate(fullMoon.getDate() - 1);

    // Find the next Sunday (Easter)
    while (fullMoon.getDay() !== 0){
        fullMoon.setDate(fullMoon.getDate() + 1);
    }

    return fullMoon;
}

// run
let now = new Date();
console.log(getEaster(now));
let jd = toJulianTime(now);

let result = getMoonPhase(jd);
console.log(result.moonIllumination);
// validation https://www.moongiant.com/phase/09/24/2020