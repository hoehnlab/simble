.. _glossary:
Glossary
=========

A list of terms used in SimBLE and their definitions.

.. glossary::

    B Cell
        A type of immune cell that produces antibodies. The protein responsible for binding in antibodies is also expressed on the surface of B cells and is called the B cell receptor (BCR). 
    
    B cell receptor
    BCR
        A protein complex found on the surface of B cells, which binds to specific antigens. The binding region of the BCR is formed by the variable regions of the heavy and light chains.

    Clone
        A group of B cells that all derive from a common ancestor B cell. In SimBLE, a clone is initiated by a single naive B cell.

    Target sequence
        The BCR amino acid sequence with theoretical "optimal" binding to an arbitrary antigen.

    Affinity
        The strength of binding between an antibody and its antigen. In SimBLE, affinity is modeled as a function of the similarity between the BCR and the target sequence.

    Germinal Center
    GC
        A specialized microenvironment within secondary lymphoid organs where B cell proliferation, somatic hypermutation, and selection occur during an immune response.

    Somatic Hypermutation
    SHM
        A cellular mechanism by which B cells mutate their BCR genes at a high rate during an immune response, leading to increased diversity in the BCR repertoire. In SimBLE, SHM is modeled using the S5F targeting model by default, but can also be run with uniform mutation in neutral simulations.

    Chain
    Heavy Chain
    Light Chain
        One of the two types of polypeptide chains that make up the BCR. In SimBLE, chain refers to both the nucleotide sequence coding for the polypeptide chain and the amino acid sequence that makes up that chain. Each BCR consists of one heavy chain (coded by the IGH locus) and one light chain (coded by either the IGL or IGK locus). 

    Migration
        The movement of B cells between different locations, such as from the germinal center to other parts of the body. In SimBLE, migration is modeled with a specified migration rate, and is closely linked to differentiation.

    Differentiation
        The process by which B cells change from one type to another, such as from germinal center B cells to memory B cells or plasma cells. In SimBLE, rates of differentiation are modeled based on time. MBCs are produced early in the germinal center reaction, while PCs are produced later. GC B cells differentiate stochastically into MBCs, while PCs arise from high-affinity MBCs. MBCs and PCs are assumed to migrate out of the germinal center upon differentiation.

    Germinal Center B Cell
    GC B Cell
        A B cell that is actively participating in the germinal center reaction, undergoing proliferation, somatic hypermutation, and selection.

    Memory B Cell
    MBC
        A long-lived B cell that has been generated during an immune response.

    Plasma Cell
    PC
        A differentiated B cell that produces and secretes large amounts of antibodies. In SimBLE, antibody production is not explicitly modeled.

    Neutral Simulation
        A simulation in which all B cells have equal fitness, meaning there is no selection based on affinity. This allows for the study of genetic drift and other neutral evolutionary processes.

    Selection
        The process by which B cells with higher affinity for the antigen are preferentially expanded and survive in the germinal center. In SimBLE, selection is modeled based on the affinity of BCRs to the target sequence.

    Uniform Neutral Simulation
        A type of neutral simulation where starting sequences are randomly generated nucleotide sequences, and all possible mutations are equally likely, without any SHM biases.