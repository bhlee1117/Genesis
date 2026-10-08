"""Genesis - a shared cloning designer for the lab's plasmid inventory.

Indexes the plasmid files on the lab share, normalises them into a searchable
inventory of plasmids and reusable parts, validates a specified construct in
silico, and generates ranked cloning routes in which every step carries an
explicit laboratory operation code.

Standard library only: no third-party packages are required to run the
importer, the design engine, the web server or the interface.

    python3 -m genesis init
    python3 -m genesis scan
    python3 -m genesis serve
"""

__version__ = "0.1.0"
__all__ = ["config", "seqio", "seqtools", "vocab", "db", "importer",
           "construct", "validate", "operations", "enzymes", "primers",
           "planner", "design", "export", "server"]
