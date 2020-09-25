/**
 * @property {number} moonPhase Phase of the moon as a fraction
 * @property {number} moonIllumination Illuminated fraction
 * @property {number} moonAgeInDays Age of moon in days
 * @property {number} distanceInKm Distance in kilometres
 * @property {number} angularDiameterInDeg Angular diameter in degrees
 * @property {number} distanceToSun Distance to sun
 * @property {number} sunAngularDiameter Sun's angular diameter
 */
export type PhaseResult = {
    moonPhase: number;
    moonIllumination: number;
    moonAgeInDays: number;
    distanceInKm: number;
    angularDiameterInDeg: number;
    distanceToSun: number;
    sunAngularDiameter: number;
}