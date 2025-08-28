class InvalidHGVSName(ValueError):
    def __init__(self, name="", part="name", reason=""):
        if name:
            message = f'Invalid HGVS {part} "{name}"'
        else:
            message = f"Invalid HGVS {part}"
        if reason:
            message += ": " + reason
        super().__init__(message)

        self.name = name
        self.part = part
        self.reason = reason


class InvalidKindError(InvalidHGVSName):
    def __init__(self, allele: str):
        super().__init__(allele, "allele", 'expected kind "c.", "p.", "g.", etc')


class NonMatchingAlleleError(InvalidHGVSName):
    def __init__(self, name: str, allele_type: str, reason="couldn't find a match"):
        super().__init__(name, allele_type + " allele", reason)


class StartGreaterThanEndError(InvalidHGVSName):
    def __init__(self, name: str):
        super().__init__(name, "coordinates", "start greater than end")
