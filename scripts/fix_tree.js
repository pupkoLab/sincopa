let removeTaxaModule = null;

async function initializeRemoveTaxa() {
    removeTaxaModule = await createRemoveTaxaModule({
        noInitialRun: true
    });

    console.log("removeTaxa WebAssembly ready");
}

initializeRemoveTaxa();


function removeBootstrapValues(treeText) {
    // Python:
    // re.sub(r'\)\d+:', '):', tree_as_str)
    return treeText.replace(/\)\d+:/g, '):');
}


function getMsaTaxa(msaText) {
    const taxa = [];

    const regex = /^>(\S+)/gm;
    let match;

    while ((match = regex.exec(msaText)) !== null) {
        taxa.push(match[1]);
    }

    return taxa;
}


/*
 * Equivalent purpose to Bio.Phylo:
 *
 * [leaf.name for leaf in tree.get_terminals()]
 *
 * For now this extracts Newick leaf names.
 */
function getTreeLabels(treeText) {

    const labels = [];

    /*
     * Leaf labels occur following "(" or ","
     * and terminate at :, comma, or ")".
     */
    const regex = /(?:\(|,)\s*([^():,;\s]+)\s*(?=[:),])/g;

    let match;

    while ((match = regex.exec(treeText)) !== null) {
        labels.push(match[1]);
    }

    return labels;
}


async function fixTree(msaText, treeText) {

    if (!removeTaxaModule) {
        throw new Error("removeTaxa WebAssembly is still loading.");
    }

    /*
     * Step 1:
     * remove bootstrap values
     */
    const treeWithoutBootstrap =
        removeBootstrapValues(treeText);


    /*
     * Step 2:
     * get taxa from tree
     */
    const treeTaxa =
        getTreeLabels(treeWithoutBootstrap);


    /*
     * Step 3:
     * get taxa from MSA
     */
    const msaTaxa =
        getMsaTaxa(msaText);

    const msaSet =
        new Set(msaTaxa);


    /*
     * Step 4:
     * determine which tree taxa should be removed
     */
    const taxaToRemove =
        treeTaxa.filter(
            taxon => !msaSet.has(taxon)
        );


    console.log("Tree taxa:", treeTaxa);
    console.log("MSA taxa:", msaTaxa);
    console.log("Taxa to remove:", taxaToRemove);


    /*
     * Step 5:
     * put files into the Emscripten filesystem
     */
    const FS = removeTaxaModule.FS;

    const inputTreePath = "/tree.no_bootstrap";
    const taxaPath = "/taxa_to_remove.txt";
    const outputTreePath = "/tree_fixed.newick";


    /*
     * Clean old files in case user runs twice.
     */
    for (const path of [
        inputTreePath,
        taxaPath,
        outputTreePath
    ]) {
        try {
            FS.unlink(path);
        }
        catch (e) {
        }
    }


    FS.writeFile(
        inputTreePath,
        treeWithoutBootstrap
    );

    FS.writeFile(
        taxaPath,
        taxaToRemove.join("\n")
    );


    /*
     * Equivalent to:
     *
     * removeTaxa input_tree taxa_file output_tree
     */
    removeTaxaModule.callMain([
        inputTreePath,
        taxaPath,
        outputTreePath
    ]);


    /*
     * Step 6:
     * retrieve fixed tree from WASM filesystem
     */
    const fixedTree =
        FS.readFile(
            outputTreePath,
            { encoding: "utf8" }
        );


    return {
        fixedTree,
        taxaToRemove,
        treeTaxa,
        msaTaxa
    };
}


document.addEventListener(
    "DOMContentLoaded",
    function () {

        const form =
            document.getElementById("form");


        form.addEventListener(
            "submit",
            async function (event) {

                event.preventDefault();

                try {

                    const msaFile =
                        document.getElementById("msa_file").files[0];

                    const treeFile =
                        document.getElementById("tree_file").files[0];


                    if (!msaFile) {
                        alert("Please select a DNA MSA file.");
                        return;
                    }

                    if (!treeFile) {
                        alert("Please select a phylogenetic tree file.");
                        return;
                    }


                    const msaText =
                        await msaFile.text();

                    const treeText =
                        await treeFile.text();


                    const result =
                        await fixTree(
                            msaText,
                            treeText
                        );


                    /*
                     * Save result for results_fix_tree.html
                     */
                    sessionStorage.setItem(
                        "fixedTree",
                        result.fixedTree
                    );

                    sessionStorage.setItem(
                        "taxaToRemove",
                        JSON.stringify(
                            result.taxaToRemove
                        )
                    );

                    sessionStorage.setItem(
                        "treeFileName",
                        treeFile.name
                    );


                    window.location.href =
                        "results_fix_tree.html";

                }
                catch (error) {

                    console.error(error);

                    alert(
                        "Tree processing failed:\n" +
                        error.message
                    );
                }
            }
        );
    }
);
