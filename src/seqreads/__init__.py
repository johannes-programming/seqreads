import dataclasses as _dataclasses

import Bio.Seq as _Seq
import Bio.SeqRecord as _SeqRecord
import na_quantors as _na


@_dataclasses.dataclass(frozen=True)
class SeqRead:
    qualities:float
    seq:_Seq.Seq
    def __post_init__(self, *args, **kwargs):
        cls = type(self)
        for k, v in __annotations__.items():
            attr = getattr(self, k)
            if type(attr) is not v:
                raise TypeError
        if len(self) and _na.notna(self.qualities):
            raise ValueError
    @classmethod
    def by_seq_and_qualities(cls, seq="", qualities='nan'):
        seq = _Seq.Seq(rec.seq)
        qualities = float(qualities)
        return cls(seq=seq, qualities=qualities)
    @classmethod
    def by_seqRecord(cls, rec):
        try:
            qualities = float(sum(rec.letter_annotations['phred_quality']))
        except:
            qualities = 'nan'
        return cls.by_seq_and_qualities(seq=rec.seq, qualities=qualities)
    @classmethod
    def by_seqReads(cls, *reads):
        ans = cls.by_seq_and_qualities("", 0)
        for read in reads:
            if type(read) is not cls:
                raise TypeError
            if len(read):
                ans += read
        return ans
    def to_record(self):
        rec = _SeqRecord.SeqRecord(self.seq)
        if len(self) == 0:
            return rec
        try:
            q = round(self.qv)
        except:
            pass
        else:
            rec.letter_annotations['phred_quality'] = [q] * len(self)
        return rec
    def __len__(self):
        return len(self.seq)
    def __add__(self, other):
        return self.by_seqReads(self, other)
    def __mul__(self, other):
        cls = type(self)
        return cls(
            self.seq * other,
            self.qualities * other,
        )
    def __rmul__(self, other):
        cls = type(self)
        return cls(
            other * self.seq,
            other * self.qualities,
        )
    @property
    def qv(self):
        try:
            return self.qualities / len(self)
        except ZeroDivisionError:
            return float('nan')