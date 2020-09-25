/**
 * @property {number} phase Phase of the moon as a fraction
 * @property {number} pphase Illuminated fraction
 * @property {number} mage Age of moon in days
 * @property {number} dist Distance in kilometres
 * @property {number} angdia Angular diameter in degrees
 * @property {number} sudist Distance to sun
 * @property {number} suangdia Sun's angular diameter
 */
export type PhaseResult = {
    phase: number;
    pphase: number;
    mage: number;
    dist: number;
    angdia: number;
    sudist: number;
    suangdia: number;
}