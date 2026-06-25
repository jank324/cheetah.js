// species.js - Named particle species

export const electron_mass_eV = 510998.95;
export const proton_mass_eV = 938272088.16;
export const deuteron_mass_eV = 1875612942.57;
export const elementary_charge = 1.602176634e-19;
export const eV_to_kg = 1.7826619216278975e-36;

export class Species {
    static known = {
        "electron": { num_elementary_charges: -1, mass_eV: electron_mass_eV },
        "positron": { num_elementary_charges: 1, mass_eV: electron_mass_eV },
        "proton": { num_elementary_charges: 1, mass_eV: proton_mass_eV },
        "antiproton": { num_elementary_charges: -1, mass_eV: proton_mass_eV },
        "deuteron": { num_elementary_charges: 1, mass_eV: deuteron_mass_eV },
    };

    constructor(name, num_elementary_charges = null, mass_eV = null) {
        this.name = name;
        if (Species.known[name]) {
            this.num_elementary_charges = Species.known[name].num_elementary_charges;
            this.mass_eV = Species.known[name].mass_eV;
        } else {
            if (num_elementary_charges === null || mass_eV === null) {
                throw new Error(`Custom species '${name}' must have charge and mass provided.`);
            }
            this.num_elementary_charges = num_elementary_charges;
            this.mass_eV = mass_eV;
        }
    }

    get mass_kg() {
        return this.mass_eV * eV_to_kg;
    }

    get charge_coulomb() {
        return this.num_elementary_charges * elementary_charge;
    }

    clone() {
        return new Species(this.name, this.num_elementary_charges, this.mass_eV);
    }
}
