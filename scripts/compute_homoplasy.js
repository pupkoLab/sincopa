let mpreconstructModule = null;


/*
 * Load MPreconstruct WebAssembly.
 */
async function initializeMPreconstruct() {
    mpreconstructModule =
        await createMPreconstructModule({
            noInitialRun: true
        });

    console.log("MPreconstruct WebAssembly ready");
}

initializeMPreconstruct();


/*
 * Browser equivalent of compute_homoplasy.py.
 *
 * Inputs:
 *   fixedMsa  - fixed MSA text
 *   fixedTree - fixed Newick tree text
 *
 * Returns:
 *   homoplasyText
 *   controlText
 */
async function computeHomoplasy(fixedMsa, fixedTree) {

    if (!mpreconstructModule) {
        throw new Error(
            "MPreconstruct WebAssembly is still loading."
        );
    }

    const FS = mpreconstructModule.FS;

    const msaPath =
        "/msa_fixed.fasta";

    const treePath =
        "/tree_fixed.newick";

    const controlPath =
        "/control.txt";

    const outputPath =
        "/homoplasy.txt";

    const logPath =
        "/control.log";

    const mpOutputPath =
        "/control.out";


    /*
     * Remove files left from a previous run.
     */
    for (const path of [
        msaPath,
        treePath,
        controlPath,
        outputPath,
        logPath,
        mpOutputPath
    ]) {
        try {
            FS.unlink(path);
        }
        catch (e) {
        }
    }


    /*
     * Put fixed inputs into the Emscripten filesystem.
     */
    FS.writeFile(
        msaPath,
        fixedMsa
    );

    FS.writeFile(
        treePath,
        fixedTree
    );


    /*
     * Equivalent to the control file generated
     * by compute_homoplasy.py.
     */
    const controlText =
`### parameters for MPreconstruct

## in and out files
_treefile ${treePath}
_seqfile ${msaPath}
_logfile ${logPath}
_outfile ${mpOutputPath}

## types are: nuc, amino, threeState, integer
_alphabetType nuc

## alphabet: a, c, g, t, and -
_alphabetSize 5

## types are: file,fitch,diff,diffSquare
_costMatrixType fitch`;


    FS.writeFile(
        controlPath,
        controlText
    );


    /*
     * Equivalent to:
     *
     * MPreconstruct control.txt homoplasy.txt
     */
    mpreconstructModule.callMain([
        controlPath,
        outputPath
    ]);


    /*
     * Read generated homoplasy result.
     */
    const homoplasyText =
        FS.readFile(
            outputPath,
            { encoding: "utf8" }
        );


    /*
     * Read the actual saved control file.
     */
    const savedControlText =
        FS.readFile(
            controlPath,
            { encoding: "utf8" }
        );


    return {
        homoplasyText,
        controlText: savedControlText
    };
}
